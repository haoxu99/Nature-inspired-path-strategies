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
#define H 0.4  //默认的路径厚度（适用于Case2）
#define With_path 0.8
#define With_material 1.75
#define pi 3.14
#define N 1	//首层与末层的层数
#define n_p 1	//频率

using namespace std;
typedef unsigned char boolean;
struct Point {
	double x, y, z, t;
	bool b;
}p;

vector<Point>coord;  //缓存一层的轮廓线

vector<Point> buffer1;
vector<vector<Point>>buffer2;
vector<vector<vector<Point>>>P;

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;//定义完整MyMesh类，用三角网格，自定的属性
MyMesh mesh;

double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;//用于生成包围盒
int X, Y, Z;//包围盒内离散点的数量

const string file_1 = "qiao.STL";

// 读取文件的函数
void readfile(string file) {//文件读取到了mesh当中
	// 请求顶点法线 vertex normals
	mesh.request_vertex_normals();
	//如果不存在顶点法线，则报错 
	if (!mesh.has_vertex_normals())
	{
		cout << "错误：标准定点属性 “法线”不存在" << endl;
		return;
	}
	// 如果有顶点发现则读取文件 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "无法读取文件:" << file << endl;
		return;
	}
	else cout << "成功读取文件:" << file << endl;
	cout << endl; // 为了ui显示好看一些
				  //如果不存在顶点法线，则计算出
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// 通过面法线计算顶点法线
		mesh.request_face_normals();
		// mesh计算出顶点法线
		mesh.update_normals();
		// 释放面法线
		mesh.release_face_normals();
	}
}

//截面函数
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm) //pt截面上一点，pnorm截面法向
{
	const float ERR = 0.001;
	//参数包括 pt，pnorm，*pilist，pnum[]  具体函数原理 见 截面算法.docx
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

	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //遍历整个模型所有的边，有交点的把id记录在flag中
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
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);//d1到面的距离
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))//线段与面相交
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];

			if (dist < mind)
			{
				nhe = hh;	//最短边
				mind = dist;//最短距离
				ne = id;    //最短所在边的编号               //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc)
	{
		delete[]flag;
		return false; //没有交点，这里return false，跳出整个函数
	}

	starte = ne;//标记循环起始的边

	suc = false;

	nhf = mesh.face_handle(nhe);//最短边所在面

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
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) //nhf最短边所在面,迭代这个面的边，只有3个
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

				MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);//获取反向半边

				nhf = mesh.face_handle(halfnext);//返回这个边所在的面

				int ne2 = halfnext.idx();

				flag[ne1] = flag[ne2] = false;//ne1,ne2是对向的两个半边，存过其中一个的交点就都变为false

				if (nhf.idx() == -1)
				{
					starte = ne;
					flag[ne] = false;
					break;
				}
				ne = ne2;//以对边的那个面再找下一个与面相交的线

				//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
				//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
				//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
				//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);

				break;
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //检索有没有第二个环
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
				p = { -10000,-10000 };//两个环中间的间隔数据
				coord.push_back(p);

				suc = false;
				break;
			}
		}


	};

	delete[]flag;
	return true;
}


//通过遍历点获取模型的包围盒
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


double distance(Point a, Point b)	//两点间距离
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2) + pow((a.z - b.z), 2));
}

//示例一路径生成
void Case1_paths() {
	Point s;
	//生成离散点，距离为路径宽度
	s.x = 0; s.y = 0; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2);

	for (int i = 0; i < Z; i++) {
		if (i == 0) {//第一层
			for (int j = 0; j < X; j++) {
				s.x += With_path;
				for (int m = 0; m < Y; m++) {
					s.y += With_path;

					(i + 1 + j + 1 + m + 1) % 2 != 0 ? s.z = MIN / 2 : s.z = MAX / 2;

					buffer1.push_back(s);
				}
				buffer2.push_back(buffer1);
				buffer1.clear();
				s.y = 0;
			}
			s.x = 0;
			P.push_back(buffer2);
			buffer2.clear();
		}
		else if (i == Z - 1) {//最后一层
			if ((i + 1) % 2 != 0) {//假如为奇数层
				for (int j = 0; j < X; j++) {
					s.x += With_path;
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(i + 1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = 0;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(i + 1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN / 2 : s.z = P[i - 1][m][j].z + MAX / 2;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = 0;
			}
		}
		else {//中间层
			if ((i + 1) % 2 != 0) {//假如为奇数层
				for (int j = 0; j < X; j++) {
					s.x += With_path;
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(i + 1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = 0;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(i + 1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MIN : s.z = P[i - 1][m][j].z + MAX;
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = 0;
			}
		}
	}

}

//示例二路径生成
void Case2_paths() {
	Point s;
	s.x = 0; s.y = 0; s.z = 0; s.t = 0;
	Z = (Bmax_z - Bmin_z) / (MAX / 2 + MIN / 2)-1;
	for (int i = 0; i < Z; i++) {
		if (i == 0) {//第一层
			for (int j = 0; j < X; j++) {
				s.x += With_path;//假如为奇数层
				for (int m = 0; m < Y; m++) {
					s.y += With_path;
					if (i == 0)
						(1 + j + 1 + m + 1) % 2 != 0 ? s.z = MIN : s.z = MAX;
					else
						(1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][j][m].z + MIN : s.z = P[i - 1][j][m].z + MAX;
					buffer1.push_back(s);
				}
				buffer2.push_back(buffer1);
				buffer1.clear();
				s.y = 0;
			}
			P.push_back(buffer2);
			buffer2.clear();
			s.x = 0;
		}
		else if (i == Z - 1) {//最后一层//假如为奇数层
			if ((i + 1) % 2 != 0) {
				for (int j = 0; j < X; j++) {
					s.x += With_path;//假如为奇数层
					for (int m = 0; m < Y; m++) {
						s.y += With_path;
						(1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//不用考虑奇偶层，应该与第一层互补
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.y = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = 0;
			}
			else {
				for (int j = 0; j < Y; j++) {
					s.y += With_path;
					for (int m = 0; m < X; m++) {
						s.x += With_path;
						(1 + j + 1 + m + 1) % 2 != 0 ? s.z = P[i - 1][m][j].z + MAX : s.z = P[i - 1][m][j].z + MIN;//不用考虑奇偶层，应该与第一层互补
						buffer1.push_back(s);
					}
					buffer2.push_back(buffer1);
					buffer1.clear();
					s.x = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = 0;
			}
		}
		else {//中间层//假如为奇数层
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
					s.y = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.x = 0;
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
					s.x = 0;
				}
				P.push_back(buffer2);
				buffer2.clear();
				s.y = 0;
			}
		}

	}
}

//判断点是否在多边形内，在true，否false，使用的是引射线法
bool InOrOutPolygon(Point a, vector<Point> polypoint) {
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// 点在两个x之间 且以点垂直y轴向上做射线
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope* (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	if (crossings % 2 != 0) {
		return true;
	}
	else {
		return false;
	}
}

//判断生成的离散点是否在模型内
void findIntersect() {
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 0; m < P[i][j].size(); m++) {
				MyMesh::Normal vf(0, 0, 1);//截面法线
				MyMesh::Point pt;//截面上一点
				pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = P[i][j][m].z;
				IntersectPlane(pt, vf);//生成一层轮廓线，存入coord
				if (!InOrOutPolygon(P[i][j][m], coord)) {
					P[i][j][m].b = false;
				}
				else {
					P[i][j][m].b = true;
				}
				coord.clear();
			}
		}
		
	}
}


void Optimization() {
	Point s;
	double k1,k2;	//斜率

	for (int i = 0; i < P.size(); i++) {//中间插入过渡插值
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {

				k1 = (P[i][j][m].z - P[i][j][m - 1].z) / With_path;
				//k2 = (P[i][j][m].t - P[i][j][m - 1].t) / With_path;
				for (int c = 1; c < With_path * 10; c++) {//插值间隔为0.1mm
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

void  Case1() {
	Case1_paths();
	findIntersect();  //前三个函数切片，产生model
	Optimization();

	double k = (MAX / 2 - MIN / 2) / With_path;//存储斜率
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				//首层
				if (i == 0) {
					P[i][j][m].t = 2 * (With_path) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(With_material, 2) * pi);
				}//末层(奇数)
				else if (i == Z - 1&&(i+1)%2!=0) {
					P[i][j][m].t = P[0][j][m].t;
				}//末层（偶数）
				else if (i == Z - 1 && (i + 1) % 2 == 0) {
					P[i][j][m].t = P[1][j][m].t / 2;
				}//中间层
				else {
					if (P[i][j][m].z > P[i][j][m - 1].z) {
						if(m % (int)(With_path * 10) == 0)
							P[i][j][m].t = 2 * (With_path) * (MAX + (MIN + 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);
						else
							P[i][j][m].t = 2 * (With_path) * ((MIN + 2 * (m % (int)(With_path * 10) * 0.1 * k)) + (MIN + 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);

					}
					else {
						if (m % (int)(With_path * 10) == 0)
							P[i][j][m].t = 2 * (With_path) * (MIN + (MAX - 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);
						else
							P[i][j][m].t = 2 * (With_path) * ((MAX - 2 * (m % (int)(With_path * 10) * 0.1 * k)) + (MAX - 2 * ((m - 1) % (int)(With_path * 10) * 0.1 * k))) / (pow(With_material, 2) * pi);
						
					}
				}
			}
		}
	}
}

void  Case2() {
	Case2_paths();
	findIntersect();  //前三个函数切片，产生model
	Optimization();

	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			for (int m = 1; m < P[i][j].size(); m++) {
				//首层
				if (i == 0) {
					P[i][j][m].t = 2 * (With_path) * (P[i][j][m].z + P[i][j][m - 1].z) / (pow(With_material, 2) * pi);
				}//末层(奇数)
				else if (i == Z - 1 && (i + 1) % 2 != 0) {
					P[i][j][m].t = 2 * (With_path) * (MAX + MIN - P[0][j][m - 1].z + MAX + MIN - P[0][j][m].z) / (pow(With_material, 2) * pi);
				} // 末层(偶数)
				else if (i == Z - 1 && (i + 1) % 2 == 0) {
					P[i][j][m].t = 2 * (With_path) * (MAX + MIN - (P[1][j][m - 1].z - H) + MAX + MIN - (P[1][j][m].z - H)) / (pow(With_material, 2) * pi);

				}//中间层
				else {
					P[i][j][m].t = 4 * (With_path) * H / (pow(With_material, 2) * pi);
				}
			}
		}
	}
}

void GcodePrint() {

	for (int i = 0; i < P.size(); i++) {//对每层中偶数条的路径翻转，为了减少跳转
		for (int j = 0; j < P[i].size(); j++) {
			if ((j + 1) % 2 == 0) {
				reverse(P[i][j].begin(), P[i][j].end());
			}
		}
	}

	FILE* fp;
	errno_t err;     //判断此文件流是否存在 存在返回1
	err = fopen_s(&fp, "ls2.gcode", "a"); //若return 1 , 则将指向这个文件的文件流给

	double E = 0;
	double r;//回抽
	int L = 50;//偏移量
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
	int o;
	for (int i = 0; i < P.size(); i++) {
		if (i > 0) {
			fprintf(fp, "G0 F3000 X%f Y%f Z%f\n", P[i][0][0].x + L, P[i][0][0].y + L, P[i][0][0].z + 1);
			fprintf(fp, "G0 F3000 X%f Y%f Z%f E%f\n", P[i][0][0].x + L, P[i][0][0].y + L, P[i][0][0].z, E += 1);
		}
		fprintf(fp, ";TYPE:FILL\n");
		for (int j = 0; j < P[i].size(); j++) {
			if (P[i][j][0].b==true)
				fprintf(fp, "G0 X%f Y%f Z%f\n", P[i][j][0].x + L, P[i][j][0].y + L, P[i][j][0].z);
			for (int m = 1; m < P[i][j].size(); m++) {
				if (P[i][j][m - 1].b==true && P[i][j][m].b==true) {
					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", P[i][j][m].x + L, P[i][j][m].y + L, P[i][j][m].z, E += distance(P[i][j][m - 1], P[i][j][m]) * P[i][j][m].t);
					o = m;
				}
				else if (P[i][j][m - 1].b == false && P[i][j][m].b == true) {
					fprintf(fp, "G0 F3000 X%f Y%f Z%f\n", P[i][j][m].x + L, P[i][j][m].y + L, P[i][j][m].z);
					fprintf(fp, "G1 F2400 X%f Y%f Z%f E%f\n", P[i][j][m].x + L, P[i][j][m].y + L, P[i][j][m].z, E += distance(P[i][j][m - 1], P[i][j][m]) * P[i][j][m].t);
				}
			}
		}
		fprintf(fp, "G0 F3000 Z%f E%f\n", P[i][P[i].size() - 1][o].z + MAX, E += -1);


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
	start = clock();//开始时间

	readfile(file_1);
	BoundingBox();//生成包围盒

	//Case1();
	Case2();

	GcodePrint();
	end = clock();//结束时间
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  //输出时间（单位：ｓ）
}
