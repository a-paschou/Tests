#include <vector>
#include <iostream>
#include "Surface.h"
#include "Line.h"
#include "Circle.h"
#include "PoissonDiskSampling.h"
#include "DelaunayTriangulation2D.h"
#include <iostream>
#include <fstream>
#include <Material.cpp>
#include <StressStrainRelation.h>
#include "CST.h"
#include <Model.h>
#include <Analyzer.h>
#include <DirectSolution.h>

using enum StressStrainRelation;
using enum DOFType;
void WriteMesh(Grid*);
void TwoTrianglesExample();
void main()
{
	double XLength{ 10 };
	double YLength{ 10 };
	double holeRadius{ 2 };
	IShape* A = new Point(0, 0);
	IShape* B = new Point(XLength, 0);
	IShape* C = new Point(XLength, YLength);
	IShape* D = new Point(0, YLength);

	IShape* AB = new Line((Point*)A, (Point*)B);
	IShape* BC = new Line((Point*)B, (Point*)C);
	IShape* CD = new Line((Point*)C, (Point*)D);
	IShape* DA = new Line((Point*)D, (Point*)A);

	Surface surface{ std::vector<IShape*> { AB, BC, CD, DA } };

	IShape* circle1 = new Circle(new Point(XLength / 2, YLength / 2), holeRadius);

	surface.Cut(std::vector<IShape*> {circle1});

	double radius{ 0.5 };
	int k{ 20 };

	PoissonDiskSampling test{ &surface, radius, k };
	IDelaunayTriangulation* mesh = new DelaunayTriangulation2D(test.PassGrid());

	WriteMesh(mesh->grid);

	// PURE TENSION PROBLEM
	double thickness{ 4 };
	Material material;
	material.elasticModulus = 21 * pow(10, 4);
	material.poissonRatio = 0.3;
	double force{ 1000 };

	// Dirichlet Boundary conditions
	std::vector<int> dirichletNodesIDs{};
	std::vector<DOFType> dirichletDofTypes{};
	std::vector<double> dirichletAmounts{};
	for (auto& node : mesh->grid->nodesMap)
	{
		if (node.second->GetX() == 0)
		{
			dirichletNodesIDs.push_back(node.second->GetID());
			dirichletDofTypes.push_back(DOFX);
			dirichletAmounts.push_back(0);
		}
	}
	for (auto& node : mesh->grid->nodesMap)
	{
		if (node.second->GetX() > 0 && node.second->GetX() < XLength && node.second->GetY() > 0 && node.second->GetY() < YLength)
		{
			dirichletNodesIDs.push_back(node.second->GetID());
			dirichletDofTypes.push_back(DOFY);
			dirichletAmounts.push_back(0);
			break;
		}
	}

	//Neumman Boundary conditions
	std::vector<int> neumannNodesIDs{};
	std::vector<DOFType> neumannDofTypes{};
	std::vector<double> neumannAmounts{};
	for (auto& node : mesh->grid->nodesMap)
	{
		if (node.second->GetX() == XLength)
		{
			neumannNodesIDs.push_back(node.second->GetID());
			neumannDofTypes.push_back(DOFX);
			neumannAmounts.push_back(0);
		}
	}

	double area{ YLength * thickness };
	double length{};
	double elementArea{};
	double elementForce{};
	for (auto& element : mesh->grid->elementsMap)
	{
		int count{};
		double y[2]{};
		int pos[2]{};
		for (size_t i = 0; i < element.second->GetNumberOfNodes(); i++)
		{
			if (mesh->grid->nodesMap.at(element.second->GetNodeID(i))->GetX() == XLength)
			{
				pos[count] = element.second->GetNodeID(i);
				y[count] = mesh->grid->nodesMap.at(element.second->GetNodeID(i))->GetY();
				count++;
			}
		}
		if (count == 2)
		{
			length = abs(y[1] - y[0]);
			elementArea = thickness * length;
			elementForce = (elementArea / area) * force;

			int nodePos{};
			for (size_t i = 0; i < neumannNodesIDs.size(); i++)
			{
				if (neumannNodesIDs.at(i) == pos[0])
				{
					neumannAmounts.at(i) += elementForce / 2;
				}
				if (neumannNodesIDs.at(i) == pos[1])
				{
					neumannAmounts.at(i) += elementForce / 2;
				}
			}
			
		}
	}

	//Build model
	Element::SetStressStrainRelations(&material, PlaneStress);
	Element::SetThickness(thickness);
	DirichletBoundaryConditions dirichlet{ &dirichletNodesIDs, &dirichletDofTypes, &dirichletAmounts };
	NeumannBoundaryConditions neumann{ &neumannNodesIDs, &neumannDofTypes, &neumannAmounts };
	BoundaryConditions boundaryConditions{ &dirichlet, &neumann };


	Model model{ mesh->grid, &boundaryConditions };
	Analyzer analyzer{ &model, new DirectSolution(FactorizationType::Cholesky) };

	int dim = model.GetLeftHandSide()->GetMatrixDimensions()[0];
	for (size_t i = 0; i < dim; i++)
	{
		double value{ model.GetLeftHandSide()->FindValue(i, i) };
		std::cout << value << std::endl;
	}

	analyzer.Solve();

	std::cout << "End" << std::endl;
}

void WriteMesh(Grid* grid)
{
	int elementType{ 2 };
	std::ofstream MyFile("mesh.msh");

	MyFile << "$MeshFormat";
	MyFile << "\n4.1 0 8 MSH4.1, ASCII";
	MyFile << "\n$EndMeshFormat";

	MyFile << "\n$Nodes";
	MyFile << "\n1 " << grid->nodesMap.size() << " 1 " << grid->nodesMap.size();
	MyFile << "\n2 1 0 " << grid->nodesMap.size();

	for (std::map<int, Node*>::iterator it = grid->nodesMap.begin(); it != grid->nodesMap.end(); ++it)
	{
		MyFile << "\n" << (it->first + 1);
	}
	for (std::map<int, Node*>::iterator it = grid->nodesMap.begin(); it != grid->nodesMap.end(); ++it)
	{
		MyFile << "\n" << it->second->GetX() << " " << it->second->GetY() << " 0";
	}
	MyFile << "\n$EndNodes";

	MyFile << "\n$Elements";
	MyFile << "\n1 " << grid->elementsMap.size() << " 1 " << grid->elementsMap.size();
	MyFile << "\n2 1 " << elementType << " " << grid->elementsMap.size();
	for (std::map<int, Element*>::iterator it = grid->elementsMap.begin(); it != grid->elementsMap.end(); it++)
	{
		MyFile << "\n" << (it->first + 1) << " " << (it->second->GetNodeID(0) + 1) << " " << (it->second->GetNodeID(1) + 1) << " " << (it->second->GetNodeID(2) + 1);
	}
	MyFile << "\n$EndElements";

	MyFile.close();
}

void TwoTrianglesExample()
{
	//Material material;
	//material.elasticModulus = 30 * pow(10, 6);
	//material.poissonRatio = 0.25;

	//Node node1{ 0,3,0,0 };
	//Node node2{ 1,3,2,0 };
	//Node node3{ 2,0,2,0 };
	//Node node4{ 3,0,0,0 };

	//std::vector<Node*> element1Nodes{ &node1, &node2, &node4 };
	//std::vector<Node*> element2Nodes{ &node3, &node4, &node2 };

	//int* nodeIDs{ nullptr };
	//nodeIDs = new int[3];
	//nodeIDs[0] = node1.GetID();
	//nodeIDs[1] = node2.GetID();
	//nodeIDs[2] = node4.GetID();
	//Element* element1 = new CST(0, nodeIDs);

	//nodeIDs[0] = node3.GetID();
	//nodeIDs[1] = node4.GetID();
	//nodeIDs[2] = node2.GetID();
	//Element* element2 = new CST(1, nodeIDs);

	//double thickness{ 0.5 };
	//Element::SetThickness(thickness);

	//Grid grid;
	//grid.nodesMap[node1.GetID()] = &node1;
	//grid.nodesMap[node2.GetID()] = &node2;
	//grid.nodesMap[node3.GetID()] = &node3;
	//grid.nodesMap[node4.GetID()] = &node4;

	//grid.elementsMap[0] = element1;
	//grid.elementsMap[1] = element2;

	//Element::SetStressStrainRelations(&material, PlaneStress);

	//std::vector<int> dirichletNodesIDs{};
	//dirichletNodesIDs.push_back(node1.GetID());
	//dirichletNodesIDs.push_back(node3.GetID());
	//dirichletNodesIDs.push_back(node3.GetID());
	//dirichletNodesIDs.push_back(node4.GetID());
	//dirichletNodesIDs.push_back(node4.GetID());

	//std::vector<DOFType> dirichletDofTypes{};
	//dirichletDofTypes.push_back(DOFY);
	//dirichletDofTypes.push_back(DOFX);
	//dirichletDofTypes.push_back(DOFY);
	//dirichletDofTypes.push_back(DOFX);
	//dirichletDofTypes.push_back(DOFY);

	//std::vector<double> dirichletAmounts{};
	//dirichletAmounts.push_back(0);
	//dirichletAmounts.push_back(0);
	//dirichletAmounts.push_back(0);
	//dirichletAmounts.push_back(0);
	//dirichletAmounts.push_back(0);
	//DirichletBoundaryConditions dirichlet{ &dirichletNodesIDs, &dirichletDofTypes, &dirichletAmounts };

	//std::vector<int> neumannNodesIDs{ node2.GetID() };
	//std::vector<DOFType> neumannDofTypes{ DOFY };
	//std::vector<double> neumannAmounts{ 1000 };
	//NeumannBoundaryConditions neumann{ &neumannNodesIDs, &neumannDofTypes, &neumannAmounts };
	//BoundaryConditions boundaryConditions{ &dirichlet, &neumann };


	//Model model{ &grid, &boundaryConditions };
	//Analyzer analyzer{ &model, new DirectSolution(FactorizationType::Cholesky) };

	//analyzer.Solve();

	//delete[] nodeIDs;
}
