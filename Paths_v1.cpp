#include<iostream>
#include<OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <string>
#include <vector>
#include <algorithm>
#include<stdio.h>
#include "clipper.hpp"  
#include"clipper.cpp"
#include<queue>
#define use_int32
using namespace ClipperLib;
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
#define MAX 0.6
#define MIN 0.2
#define H 0.4  //Ĭ�ϵ�·����ȣ�������Case2��
#define With_path 0.8*2
#define With_material 1.75
#define pi 3.14
#define N 1	//�ײ���ĩ��Ĳ���
#define NP 4	//��������

#define V_size 2 //����֧���㷨�����ش�С
#define With_support 0.4 //֧�ŵ�·�����

using namespace std;
typedef unsigned char boolean;
struct Point {
	double x, y, z, t;
	bool b = false;
};

vector<Point>coord;  //����һ���������

vector<Point> buffer1;
vector<vector<Point>>buffer2;
vector<vector<vector<Point>>>paths;
vector<vector<vector<Point>>>P;

vector<vector<vector<Point>>>P_v;//֧���㷨�����ص�
vector<vector<vector<Point>>>Path_support;//֧���㷨��·��

vector<vector<Point>>model; //slice  ����ģ�͵ģ��ֲ�����������꣬һά��������ά�˲�������ߵ�����
//vector<vector<Point>>model2; //slice  ����ģ�͵ģ��ֲ�����������꣬һά��������ά�˲�������ߵ�����
vector<Paths>model3;//����Ƭ��model��ת�����͵õ��ģ���ƫ��֮����������
vector<Point>shellpath;//������ϵĵ�
vector<vector<Point> >shellpaths;//��ĳһ������ϵ�·��
vector<vector<vector<Point> > >shell1;//������������ϵ�·��      һά��������άÿ�㲻ͬ�Ľ���·������ά·���ĵ�����

vector<Paths>model2;//����Ƭ��model��ת�����͵õ��ģ���ƫ��֮����������
vector<vector<vector<Point> > >shell2;     //��ĳһ������ϵ�·������  һά�ò�仯����������ά�ǵĲ���·������ά·������


vector<vector<vector<vector<Point>>>> shell;//����������ϵ�·��    һά��������άÿ�㲻ͬ�Ľ���·������ά��ͬ����·���ǵĲ�������ά��·���������


typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;//��������MyMesh�࣬�����������Զ�������
MyMesh mesh;

double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;//�������ɰ�Χ��
int X, Y, Z;//��Χ������ɢ�������

const string file_1 = "Fang.STL";

// ��ȡ�ļ��ĺ���
void readfile(string file) {//�ļ���ȡ����mesh����
	// ���󶥵㷨�� vertex normals
	mesh.request_vertex_normals();
	//��������ڶ��㷨�ߣ��򱨴� 
	if (!mesh.has_vertex_normals())
	{
		cout << "���󣺱�׼�������� �����ߡ�������" << endl;
		return;
	}
	// ����ж��㷢�����ȡ�ļ� 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "�޷���ȡ�ļ�:" << file << endl;
		return;
	}
	else cout << "�ɹ���ȡ�ļ�:" << file << endl;
	cout << endl; // Ϊ��ui��ʾ�ÿ�һЩ
				  //��������ڶ��㷨�ߣ�������
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// ͨ���淨�߼��㶥�㷨��
		mesh.request_face_normals();
		// mesh��������㷨��
		mesh.update_normals();
		// �ͷ��淨��
		mesh.release_face_normals();
	}
}

//���溯��
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm) //pt������һ�㣬pnorm���淨��
{
	const float ERR = 0.001;
	//�������� pt��pnorm��*pilist��pnum[]  ���庯��ԭ�� �� �����㷨.docx
	int starte, ne, ne1, nf;
	MyMesh::Point vt1, vt2;
	//MyMesh::Face f1;
	MyMesh::HalfedgeHandle nhe;
	MyMesh::FaceHandle nhf;
	float d1, d2, sd1, sd2;
	bool* flag, suc;
	float dist, mind = 1.0e+8;

	sd1 = sd2 = -10000;
	int esize = mesh.n_halfedges();
	flag = new bool[esize];
	Point p;
	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //��������ģ�����еıߣ��н���İ�id��¼��flag��
	{

		MyMesh::HalfedgeHandle hh = *it;
		int id = hh.idx();
		flag[id] = false;

		auto fromVertex = mesh.from_vertex_handle(hh);
		auto toVertex = mesh.to_vertex_handle(hh);
		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);
		//printf("$ %.3f %.3f $\n", vt1.data()[0],vt2.data()[0]);
		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);//d1����ľ���
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))//�߶������ཻ
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];

			if (dist < mind)
			{
				nhe = hh;	//��̱�
				mind = dist;//��̾���
				ne = id;    //������ڱߵı��               //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc)
	{
		delete[]flag;
		return false; //û�н��㣬����return false��������������
	}

	starte = ne;//���ѭ����ʼ�ı�

	suc = false;

	nhf = mesh.face_handle(nhe);//��̱�������

	while (!suc)
	{
		//printf("%%%%");	

		auto fromVertex = mesh.from_vertex_handle(nhe);
		auto toVertex = mesh.to_vertex_handle(nhe);

		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);

		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);
		//printf("$$$%lf %lf \n", d1, d2);
		if ((sd1 == d1) && (sd2 == d2))
		{
			flag[ne] = false;
		}
		sd1 = d1; sd2 = d2;
		//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
		//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
		//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];

		p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
		coord.push_back(p);


		do {
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) //nhf��̱�������,���������ıߣ�ֻ��3��
			{
				MyMesh::HalfedgeHandle halfnow = *it;

				const int ne1 = halfnow.idx();

				if (flag[ne1] == false || ne == ne1) continue;

				MyMesh::VertexHandle fromV = mesh.from_vertex_handle(halfnow);
				MyMesh::VertexHandle toV = mesh.to_vertex_handle(halfnow);

				vt1 = mesh.point(fromV);
				vt2 = mesh.point(toV);

				d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
				d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

				p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2)) * (vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
				coord.push_back(p);

				MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);//��ȡ������

				nhf = mesh.face_handle(halfnext);//������������ڵ���

				int ne2 = halfnext.idx();

				flag[ne1] = flag[ne2] = false;//ne1,ne2�Ƕ����������ߣ��������һ���Ľ���Ͷ���Ϊfalse

				if (nhf.idx() == -1)
				{
					starte = ne;
					flag[ne] = false;
					break;
				}
				ne = ne2;//�ԶԱߵ��Ǹ���������һ�������ཻ����

				//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
				//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
				//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
				//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);

				break;
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //������û�еڶ�����
		{

			MyMesh::HalfedgeHandle hh = *it;

			int id = hh.idx();

			if (flag[id] == true)
			{
				ne = id;
				starte = ne;
				nhe = hh;
				nhf = mesh.face_handle(nhe);
				if (nhf.idx() == -1)
				{
					flag[ne] = false;
					continue;
				}
				//pilist[num].data()[0] = -10000;
				//pilist[num].data()[1] = -10000;
				//pilist[num].data()[2] = -10000;
				p = { -10000,-10000 };//�������м�ļ������
				coord.push_back(p);

				suc = false;
				break;
			}
		}


	};

	delete[]flag;
	return true;
}


//ͨ���������ȡģ�͵İ�Χ��
void BoundingBox() {
	MyMesh::Point pt;
	int st = 0;
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		MyMesh::VertexHandle vh_i = *it;
		pt = mesh.point(vh_i);
		px = pt.data()[0];
		py = pt.data()[1];
		pz = pt.data()[2];
		if (st == 0) {
			Bmax_x = Bmin_x = px;
			Bmax_y = Bmin_y = py;
			Bmax_z = Bmin_z = pz;
			st++;
		}
		else {
			if (px > Bmax_x)Bmax_x = px; else if (px < Bmin_x)Bmin_x = px;
			if (py > Bmax_y)Bmax_y = py; else if (py < Bmin_y)Bmin_y = py;
			if (pz > Bmax_z)Bmax_z = pz; else if (pz < Bmin_z)Bmin_z = pz;
		}
	}
	Bmin_x -= With_path / 2;
	Bmin_y -= With_path / 2;
	X = (Bmax_x - Bmin_x) / With_path;
	Y = (Bmax_y - Bmin_y) / With_path;

}


double distance(Point a, Point b)	//��������
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2) + pow((a.z - b.z), 2));
}

//ʾ��һ·������
void Case1_paths() {
	Point s;
	//������ɢ�㣬����Ϊ·�����
	s.x = Bmin_x; s.y = Bmin_y; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2) - 1;

	for (int i = 0; i < Z; i++) {
		if (i == 0) {//��һ��
			for (int j = 0; j < X; j++) {
				s.x += With_path;
				for (int m = 0; m < Y; m++) {
					s.y += With_path;

					(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = Bmin_z + MIN / 2 : s.z = Bmin_z + MAX / 2;

					buffer1.push_back(s);
				}
				buffer2.push_back(buffer1);
				buffer1.clear();
				s.y = Bmin_y;
			}
			s.x = Bmin_x;
			P.push_back(buffer2);
			buffer2.clear();
		}
		else if (i == Z - 1) {//���һ��
			if ((i + 1) % 2 != 0) {//����Ϊ������
				for (int j = 0; j < X; j++) {
					s.x += With_path;
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else {//�м��
			if ((i + 1) % 2 != 0) {//����Ϊ������
				for (int j = 0; j < X; j++) {
					s.x += With_path;
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(i + 1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
	}

}

//ʾ����·������
void Case2_paths() {
	Point s;
	s.x = Bmin_x; s.y = Bmin_y; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2) - 1;
	for (int i = 0; i < Z; i++) {
		if (i < N) {//�ײ�
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += With_path;//����Ϊ������
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = Bmin_z + MIN : s.z = Bmin_z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;//����Ϊż����
					for (int m = 0; m < X; m++) {
						s.x += With_path;
					
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else if (i >= Z - N) {//ĩ��//����Ϊ������
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += With_path;//����Ϊ������
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//���ÿ�����ż�㣬Ӧ�����һ�㻥��
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(1 + j / NP + 1 + m / NP + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//���ÿ�����ż�㣬Ӧ�����һ�㻥��
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}
		else {//�м��//����Ϊ������
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += With_path;
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						s.z = P[i - 1][m][j].z + H;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = Bmin_y;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = Bmin_x;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						s.z = P[i - 1][m][j].z + H;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = Bmin_x;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = Bmin_y;
			}
		}

	}
}

//�жϵ��Ƿ��ڶ�����ڣ���true����false��ʹ�õ��������߷�
bool InOrOutPolygon(Point a, vector<Point> polypoint) {
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// ��������x֮�� ���Ե㴹ֱy������������
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope* (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	if (crossings % 2 != 0&&crossings!=0) {
		return true;
	}
	else {
		return false;
	}
}

//�ж����ɵ���ɢ���Ƿ���ģ����
void findIntersect() {
	Path  way;
	Paths layer;//ÿ��·������
	IntPoint p1;
	Point p;
	for (int i = 0; i < P.size(); i++) {
		printf("���㷨��%.2lf%%\r", i * 100.0 / P.size());
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 0; m < P[i][j].size(); m++) {
				MyMesh::Normal vf(0, 0, 1);//���淨��
				MyMesh::Point pt;//������һ��
				pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P[i][j][m].z;
				IntersectPlane(pt, vf);//����һ�������ߣ�����coord

				for (int k = 0; k < coord.size(); k++)
				{
					if (coord[k].x != -10000)
					{
						p1.X = coord[k].x * 100;
						p1.Y = coord[k].y * 100;
						way.push_back(p1);
					}
					else
					{
						layer.push_back(way);
						way.clear();
					}
				}
				layer.push_back(way);
				way.clear();
				coord.clear();

				for (int k1 = 0; k1 < layer.size(); k1++)
				{
					for (int k = 0; k < layer[k1].size(); k++)
					{
						p.x = double(layer[k1][k].X) / 100;
						p.y = double(layer[k1][k].Y) / 100;

						shellpath.push_back(p);
					}
					shellpaths.push_back(shellpath);
					shellpath.clear();
				}
				layer.clear();

				for (int k = 0; k < shellpaths.size(); k++) {
					if (InOrOutPolygon(P[i][j][m], shellpaths[k])) {
						P[i][j][m].b = true;
					}
				}
				shellpaths.clear();
				
			}
		}
		
	}
}

void Contour() {
	double z = Bmin_z;
	for (int i = 0; i < Z; i++) {
		z += (MIN + MAX) / 2;
		MyMesh::Normal vf(0, 0, 1);//���淨��
		MyMesh::Point pt; //������һ��
		pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = z;
		IntersectPlane(pt, vf);//����һ�������ߣ�����coord
		model.push_back(coord);
		coord.clear();
	}
	Path  way;
	Paths layer;//ÿ��·������
	IntPoint p1;
	ClipperOffset co;
	Paths solution;
	Point p;
	for (int i = 0; i < model.size(); i++)//����model����i��
	{
	

		for (int j = 0; j < model[i].size(); j++)//����model����i���еĵ�j����
		{
			if (model[i][j].x != -10000)
			{
				p1.X = model[i][j].x * 100;
				p1.Y = model[i][j].y * 100;
				way.push_back(p1);
			}
			else
			{
				layer.push_back(way);
				way.clear();
			}
		}
		layer.push_back(way);
		way.clear();

		model2.push_back(layer);
		layer.clear();
//*********************������Ȧ�ĵ㸳��shell1***********************************
		for (int j = 0; j < model2[i].size(); j++)
		{
			for (int k = 0; k < model2[i][j].size(); k++)
			{
				p.x = double(model2[i][j][k].X) / 100;
				p.y = double(model2[i][j][k].Y) / 100;

				shellpath.push_back(p);
			}
			shellpaths.push_back(shellpath);
			shellpath.clear();
		}
		shell1.push_back(shellpaths);
		shellpaths.clear();



		//******************************��ƫ��***********************************************
		printf("%.2lf%%\r", i * 100.0 / model.size());
		for (int j = 0; j < model2[i].size(); j++)//��ʱ��������model2��ĳ����ĵ�j��path��ƫ��
		{
			double area = abs(Area(model2[i][j]) / 10000);
			co.Clear();
			co.AddPath(model2[i][j], jtRound, etClosedPolygon);   //����׼��ƫ�Ƶ�·��
			int times = 1;	//ƫ�ƴ���������Ǵ�ӡ�Ĳ���
			int sh = With_path * 100;     //ƫ�ƺ��

			for (int a = 0; a < times; a++)					 //����ƫ��
			{

				co.Execute(solution, sh);							//ÿ��ƫ�Ƶĺ��

				for (int n = 0; n < solution[0].size(); n++)
				{
					p.x = double(solution[0][n].X) / 100;
					p.y = double(solution[0][n].Y) / 100;

					shellpath.push_back(p);
				}

				p.x = double(solution[0][0].X) / 100;
				p.y = double(solution[0][0].Y) / 100;

				shellpath.push_back(p);

				co.Clear();
				co.AddPath(solution[0], jtRound, etClosedPolygon);  //׼����һ��ƫ��
				solution.clear();

				shellpaths.push_back(shellpath);
				shellpath.clear();
			}

			shell2.push_back(shellpaths);
			shellpaths.clear();

		}
		//solution.clear();
		shell.push_back(shell2);
		shell2.clear();
	}

}

void Optimization() {
	Point s;
	double k1,k2;	//б��

	for (int i = 0; i < P.size(); i++) {//�м������ɲ�ֵ
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				if (P[i][j][m].z != P[i][j][m - 1].z) {
					k1 = (P[i][j][m].z - P[i][j][m - 1].z) / With_path;
					//k2 = (P[i][j][m].t - P[i][j][m - 1].t) / With_path;
					for (int c = 1; c < With_path * 10; c++) {//��ֵ���Ϊ0.1mm
						P[i][j][m - 1].x == P[i][j][m].x ? s.y = P[i][j][m - 1].y + 0.1 : s.x = P[i][j][m - 1].x + 0.1;
						P[i][j][m - 1].x == P[i][j][m].x ? s.x = P[i][j][m - 1].x : s.y = P[i][j][m - 1].y;
						s.z = P[i][j][m - 1].z + 0.1 * k1;
						//s.t = P[i][j][m - 1].t + 0.1 * k2;
						(P[i][j][m - 1].b == true && P[i][j][m].b == true) ? s.b = true : s.b = false;
						P[i][j].insert(P[i][j].begin() + m, s);
						m += 1;
					}
				}

			}
		}
	}

}

void  Case1() {
	Case1_paths();
	findIntersect();  //ǰ����������Ƭ������model
	Optimization();

	double k = (MAX / 2 - MIN / 2) / With_path;//�洢б��
	double z = (MIN + MAX) / 2;
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {

				//�ײ�
				if (i == 0) {
					P[i][j][m].t = 2 * (With_path) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(With_material, 2) * pi);
				}
				else if (i == 1) {//�ڶ���
					if (P[i][j][m].z > P[i][j][m - 1].z) {
						if (m % (int)(With_path * 10) == 0)
							P[i][j][m].t = 2 * (With_path) * (MAX + (MIN + 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);
						else
							P[i][j][m].t = 2 * (With_path) * ((MIN + 2 * (m % (int)(With_path * 10) * 0.1 * k)) + (MIN + 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);

					}
					else if (P[i][j][m].z < P[i][j][m - 1].z) {
						if (m % (int)(With_path * 10) == 0)
							P[i][j][m].t = 2 * (With_path) * (MIN + (MAX - 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);
						else
							P[i][j][m].t = 2 * (With_path) * ((MAX - 2 * (m % (int)(With_path * 10) * 0.1 * k)) + (MAX - 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);

					}
					else {
						if (P[i][j][m].z == MIN / 2 + MAX)
							P[i][j][m].t = 4 * (With_path)*MAX / (pow(With_material, 2) * pi);
						else
							P[i][j][m].t = 4 * (With_path)*MIN / (pow(With_material, 2) * pi);
					}
				}//ĩ��(����)
				else if (i == Z - 1 && (i + 1) % 2 != 0) {
					P[i][j][m].t = P[0][j][m].t;
				}//ĩ�㣨ż����
				else if (i == Z - 1 && (i + 1) % 2 == 0) {
					P[i][j][m].t = P[1][j][m].t / 2;
				}//�м��
				else {
					if ((i + 1) % 2 != 0) {
						P[i][j][m].t = P[0][j][m].t * 2;
					}
					else {
						P[i][j][m].t = P[1][j][m].t;
					}

				}

			}
		}
	}
}

void  Case2() {
	Case2_paths();
	findIntersect();  //ǰ����������Ƭ������model
	Optimization();

	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				//�ײ�
				if (i < N) {
					if (i == 0) {
						P[i][j][m].t = 2 * (With_path) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(With_material, 2) * pi);
					}
					else {
						P[i][j][m].t = 2 * (With_path) * (P[i][j][m].z/2 + P[i][j][m - 1].z/2) / (pow(With_material, 2) * pi);
					}
				}
				//ĩ��(����)
				else if (i >= Z - N && (i + 1) % 2 != 0) {
					P[i][j][m].t = 2 * (With_path) * (MAX + MIN - P[0][j][m - 1].z + MAX + MIN - P[0][j][m].z) / (pow(With_material, 2) * pi);
				} // ĩ��(ż��)
				else if (i >= Z - N && (i + 1) % 2 == 0) {
					P[i][j][m].t = 2 * (With_path) * (MAX + MIN - (P[1][j][m - 1].z - H) + MAX + MIN - (P[1][j][m].z - H)) / (pow(With_material, 2) * pi);

				}//�м��
				else {
					P[i][j][m].t = 4 * (With_path) * H / (pow(With_material, 2) * pi);
				}
			}
		}
	}
}

void Paths_Optimization() {
	Point s;
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 0; m < P[i][j].size(); m++) {
				if (P[i][j][m].b == true) {
					buffer1.push_back(P[i][j][m]);
				}
				else if (P[i][j][m].b == false && buffer1.size() != 0) {
					buffer2.push_back(buffer1);
					buffer1.clear();
				}
			}
			if (buffer1.size() != 0) {
				buffer2.push_back(buffer1);
				buffer1.clear();
			}
		}
		paths.push_back(buffer2);
		buffer2.clear();
	}

	double arcs, s_begin, s_end;
	for (int i = 0; i < paths.size(); i++) {
		for (int j = 0; j < paths[i].size() - 1; j++) {
			arcs = distance(paths[i][j][paths[i][j].size() - 1], paths[i][j + 1][0]);
			for (int m = j + 1; m < paths[i].size(); m++) {
				s_begin = distance(paths[i][j][paths[i][j].size() - 1], paths[i][m][0]);
				s_end = distance(paths[i][j][paths[i][j].size() - 1], paths[i][m][paths[i][m].size() - 1]);
				if (s_begin < arcs && s_begin < s_end) {
					arcs = s_begin;
					buffer1 = paths[i][m];
					paths[i][m] = paths[i][j + 1];
					paths[i][j + 1] = buffer1;
					buffer1.clear();
				}
				else if (s_end < arcs && s_end < s_begin) {
					arcs = s_end;
					reverse(paths[i][m].begin(), paths[i][m].end());
					buffer1 = paths[i][m];
					paths[i][m] = paths[i][j + 1];
					paths[i][j + 1] = buffer1;
					buffer1.clear();
				}
			}
		}
	}

}

void ContourGcodePrint() {


	FILE* fp;
	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1
	err = fopen_s(&fp, "ls2.gcode", "a"); //��return 1 , ��ָ������ļ����ļ�����
	double t = (4 * ((MIN + MAX) / 2) * With_path) / (pow(With_material, 2) * pi);
	double E = 0;
	double r;//�س�
	int L = 50;//ƫ����
	fprintf(fp, ";FLAVOR:Marlin\n");
	fprintf(fp, ";Generated with Cura_SteamEngine 4.10.0\n");
	fprintf(fp, "M140 S50\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M190 S50\n");
	fprintf(fp, "M104 S210\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M201 X500.00 Y500.00 Z100.00 E5000.00 ;Setup machine max acceleration\n");
	fprintf(fp, "M203 X500.00 Y500.00 Z10.00 E50.00 ;Setup machine max feedrate\n");
	fprintf(fp, "M204 P500.00 R1000.00 T500.00 ;Setup Print/Retract/Travel acceleration\n");
	fprintf(fp, "M205 X8.00 Y8.00 Z0.40 E5.00 ;Setup Jerk\n");
	fprintf(fp, "M220 S100 ;Reset Feedrate\n");
	fprintf(fp, "M221 S100 ;Reset Flowrate\n");

	fprintf(fp, "G28 ;Home\n");

	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F2700 E-5\n");
	fprintf(fp, "M107\n");

	double z = Bmin_z;

	for (int i = 0; i < paths.size(); i++) {
		fprintf(fp, ";LAYER:%d\n", i);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "M73 P%.f\n", float((100 * i) / paths.size()));
		z += (MIN + MAX) / 2;
		fprintf(fp, "G0 X%f Y%f Z%f\n", shell1[i][0][0].x + L, shell1[i][0][0].y + L, z);
		if (i > 0) {
			fprintf(fp, "G0 E%f\n", E += 1);
		}
		for (int j = 0; j < shell1[i].size(); j++) {
			fprintf(fp, "G0 X%f Y%f Z%f\n", shell1[i][j][0].x + L, shell1[i][j][0].y + L, z);
			 
			fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", shell1[i][j][1].x + L, shell1[i][j][1].y + L,z, E += distance(shell1[i][j][0], shell1[i][j][1]) * t);
			for (int k = 2; k < shell1[i][j].size(); k++)

			{
				fprintf(fp, "G1 X%.3f Y%.3f Z%f E%.5f\n", shell1[i][j][k].x + L, shell1[i][j][k].y + L,z, E += distance(shell1[i][j][k - 1], shell1[i][j][k]) * t);

			}
			for (int k = 0; k < shell[i][j].size(); k++)
			{
				fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.1f\n", shell[i][j][k][0].x + L, shell[i][j][k][0].y + L, z);

				fprintf(fp, "G1 F2400 X%.3f Y%.3f Z%f E%.5f\n", shell[i][j][k][1].x + L, shell[i][j][k][1].y + L,z, E += distance(shell[i][j][k][0], shell[i][j][k][1]) * t);

				for (int m = 2; m < shell[i][j][k].size(); m++)

				{
					fprintf(fp, "G1 X%.3f Y%.3f Z%f E%.5f\n", shell[i][j][k][m].x + L, shell[i][j][k][m].y + L,z, E += distance(shell[i][j][k][m - 1], shell[i][j][k][m]) * t);

				}
			}
		}
		
		fprintf(fp, ";TYPE:FILL\n");
		for (int j = 0; j < paths[i].size(); j++) {
				fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][j][0].x + L, paths[i][j][0].y + L, paths[i][j][0].z);
			for (int m = 1; m < paths[i][j].size(); m++) {
				
					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, paths[i][j][m].z, E += distance(paths[i][j][m - 1], paths[i][j][m]) * paths[i][j][m].t);
					
			}
		}
		fprintf(fp, "G0 F3000 Z%f E%f\n", paths[i][paths[i].size() - 1][paths[i][paths[i].size() - 1].size() - 1].z + MAX, E += -1);


	}
	fprintf(fp, "M140 S0\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G1 E-2 F2700\n");
	fprintf(fp, "G1 E-2 Z0.2 F2400 ;Retract and raise Z\n");
	fprintf(fp, "G1 X5 Y5 F3000 ;Wipe out\n");
	fprintf(fp, "G1 Z10 ;Raise Z more\n");
	fprintf(fp, "G90 ;Absolute positioning\n");

	fprintf(fp, "G1 X0 Y300 ;Present print\n");
	fprintf(fp, "M106 S0 ;Turn-off fan\n");
	fprintf(fp, "M104 S0 ;Turn-off hotend\n");
	fprintf(fp, "M140 S0 ;Turn-off bed\n");

	fprintf(fp, "M84 X Y E ;Disable all steppers but Z\n");

	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fclose(fp);
}

void Support() {

	Point p;
	//������Χ�����ؿ�
	for (double i = Bmin_z; i < Bmax_z; i += H) {
		for (double j = Bmin_y + V_size / 2; j < Bmax_y + V_size / 2; j += V_size) {
			for (double m = Bmin_x + V_size/2; m < Bmax_x + V_size/2; m += V_size) {
				p.x = m; p.y = j; p.z = i;
				p.t = 2 * With_support * V_size / (pow(With_material, 2) * pi);
				p.b = true;
				buffer1.push_back(p);
			}
			buffer2.push_back(buffer1);
			buffer1.clear();
		}
		P_v.push_back(buffer2);
		buffer2.clear();
	}

	//�ж�������ĵ�
	Path  way;
	Paths layer;//ÿ��·������
	IntPoint p1;
	for (int i = 0; i < P_v.size(); i++) {
		printf("���㷨��%.2lf%%\r", i * 100.0 / P_v.size());
		for (int j = 0; j < P_v[i].size(); j++) {
			for (int m = 0; m < P_v[i][j].size(); m++) {
				MyMesh::Normal vf(0, 0, 1);//���淨��
				MyMesh::Point pt;//������һ��
				pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P_v[i][j][m].z;
				IntersectPlane(pt, vf);//����һ�������ߣ�����coord

				for (int k = 0; k < coord.size(); k++)
				{
					if (coord[k].x != -10000)
					{
						p1.X = coord[k].x * 100;
						p1.Y = coord[k].y * 100;
						way.push_back(p1);
					}
					else
					{
						layer.push_back(way);
						way.clear();
					}
				}
				layer.push_back(way);
				way.clear();
				coord.clear();

				for (int k1 = 0; k1 < layer.size(); k1++)
				{
					for (int k = 0; k < layer[k1].size(); k++)
					{
						p.x = double(layer[k1][k].X) / 100;
						p.y = double(layer[k1][k].Y) / 100;

						shellpath.push_back(p);
					}
					shellpaths.push_back(shellpath);
					shellpath.clear();
				}
				layer.clear();

				for (int k = 0; k < shellpaths.size(); k++) {
					if (InOrOutPolygon(P_v[i][j][m], shellpaths[k])) {
						P_v[i][j][m].b = false;
					}
				}
				shellpaths.clear();

			}
		}

	}
	//����֧��·��
	Point s;
	for (int i = 0; i < P_v.size(); i++) {
		for (int j = 0; j < P_v[i].size(); j++) {
			for (int m = 0; m < P_v[i][j].size(); m++) {
				if (P_v[i][j][m].b == true) {
					buffer1.push_back(P_v[i][j][m]);
				}
				else if (P_v[i][j][m].b == false && buffer1.size() != 0) {
					buffer2.push_back(buffer1);
					buffer1.clear();
				}
			}
			if (buffer1.size() != 0) {
				buffer2.push_back(buffer1);
				buffer1.clear();
			}
		}
		Path_support.push_back(buffer2);
		buffer2.clear();
	}

	double arcs, s_begin, s_end;
	for (int i = 0; i < Path_support.size(); i++) {
		for (int j = 0; j < Path_support[i].size() - 1; j++) {
			arcs = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][j + 1][0]);
			for (int m = j + 1; m < Path_support[i].size(); m++) {
				s_begin = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][m][0]);
				s_end = distance(Path_support[i][j][Path_support[i][j].size() - 1], Path_support[i][m][Path_support[i][m].size() - 1]);
				if (s_begin < arcs && s_begin < s_end) {
					arcs = s_begin;
					buffer1 = Path_support[i][m];
					Path_support[i][m] = Path_support[i][j + 1];
					Path_support[i][j + 1] = buffer1;
					buffer1.clear();
				}
				else if (s_end < arcs && s_end < s_begin) {
					arcs = s_end;
					reverse(Path_support[i][m].begin(), Path_support[i][m].end());
					buffer1 = Path_support[i][m];
					Path_support[i][m] = Path_support[i][j + 1];
					Path_support[i][j + 1] = buffer1;
					buffer1.clear();
				}
			}
		}
	}


}

void GcodePrint() {
	FILE* fp;
	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1
	err = fopen_s(&fp, "ls2.gcode", "a"); //��return 1 , ��ָ������ļ����ļ�����
	double t = (4 * ((MIN + MAX) / 2) * With_path) / (pow(With_material, 2) * pi);
	double E = 0;
	double r;//�س�
	int L = 50;//ƫ����
	fprintf(fp, ";FLAVOR:Marlin\n");
	fprintf(fp, ";Generated with Cura_SteamEngine 4.10.0\n");
	fprintf(fp, "M140 S50\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M190 S50\n");
	fprintf(fp, "M104 S210\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M201 X500.00 Y500.00 Z100.00 E5000.00 ;Setup machine max acceleration\n");
	fprintf(fp, "M203 X500.00 Y500.00 Z10.00 E50.00 ;Setup machine max feedrate\n");
	fprintf(fp, "M204 P500.00 R1000.00 T500.00 ;Setup Print/Retract/Travel acceleration\n");
	fprintf(fp, "M205 X8.00 Y8.00 Z0.40 E5.00 ;Setup Jerk\n");
	fprintf(fp, "M220 S100 ;Reset Feedrate\n");
	fprintf(fp, "M221 S100 ;Reset Flowrate\n");

	fprintf(fp, "G28 ;Home\n");

	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F2700 E-5\n");
	fprintf(fp, "M107\n");

	double z = Bmin_z;

	for (int i = 0; i < paths.size(); i++) {
		fprintf(fp, ";LAYER:%d\n", i);
	
		fprintf(fp, "M73 P%.f\n", float((100 * i) / paths.size()));
		fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][0][0].x + L, paths[i][0][0].y + L, paths[i][0][0].z);
		fprintf(fp, "G0 E%f\n", E += 1);

		fprintf(fp, ";TYPE:FILL\n");
		for (int j = 0; j < paths[i].size(); j++) {
			
			fprintf(fp, "G0 X%f Y%f Z%f\n", paths[i][j][0].x + L, paths[i][j][0].y + L, paths[i][j][0].z);
			if (i > 0) {
				
			}
			for (int m = 1; m < paths[i][j].size(); m++) {

				fprintf(fp, "G1 F1200 X%f Y%f Z%f E%f\n", paths[i][j][m].x + L, paths[i][j][m].y + L, paths[i][j][m].z, E += distance(paths[i][j][m - 1], paths[i][j][m]) * paths[i][j][m].t);

			}

		}

		//support
		/*if (i > 0) {
			for (int j = 0; j < Path_support[i].size(); j++) {

				fprintf(fp, "G0 X%f Y%f Z%f\n", Path_support[i][j][0].x + L, Path_support[i][j][0].y + L, Path_support[i][j][0].z);
				if (i > 0) {

				}
				for (int m = 1; m < Path_support[i][j].size(); m++) {

					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", Path_support[i][j][m].x + L, Path_support[i][j][m].y + L, Path_support[i][j][m].z, E += distance(Path_support[i][j][m - 1], Path_support[i][j][m]) * Path_support[i][j][m].t);

				}

			}
		}*/

		fprintf(fp, "G0 F3000 Z%f E%f\n", paths[i][paths[i].size() - 1][paths[i][paths[i].size() - 1].size() - 1].z + MAX, E -= 1);


	}
	fprintf(fp, "M140 S0\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G1 E-2 F2700\n");
	fprintf(fp, "G1 E-2 Z0.2 F2400 ;Retract and raise Z\n");
	fprintf(fp, "G1 X5 Y5 F3000 ;Wipe out\n");
	fprintf(fp, "G1 Z10 ;Raise Z more\n");
	fprintf(fp, "G90 ;Absolute positioning\n");

	fprintf(fp, "G1 X0 Y300 ;Present print\n");
	fprintf(fp, "M106 S0 ;Turn-off fan\n");
	fprintf(fp, "M104 S0 ;Turn-off hotend\n");
	fprintf(fp, "M140 S0 ;Turn-off bed\n");

	fprintf(fp, "M84 X Y E ;Disable all steppers but Z\n");

	fprintf(fp, "M82 ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fclose(fp);
}

void  main(int argc, char** argv) {
	clock_t start, end;
	start = clock();//��ʼʱ��

	readfile(file_1);
	BoundingBox();//���ɰ�Χ��

	Case1();
	//Case2();
	Paths_Optimization();
	Contour();
	//Support();
	GcodePrint();
	//ContourGcodePrint();
	end = clock();//����ʱ��
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  //���ʱ�䣨��λ����
}
