#include <igl/fast_winding_number.h>
#include <igl/read_triangle_mesh.h>
#include <igl/slice_mask.h>
#include <Eigen/Geometry>
#include <igl/octree.h>
#include <igl/barycenter.h>
#include <igl/knn.h>
#include <igl/voxel_grid.h>
#include <igl/marching_cubes.h>
#include <igl/random_points_on_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/per_vertex_normals.h>
#include <igl/copyleft/cgal/point_areas.h>
#include <igl/hausdorff.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
using namespace std;
using namespace Eigen;
double interpolant(Vector3d b, VectorXd x, MatrixXd P){
  int prows=P.rows();
  double p=x(prows)+x(prows+1)*b(0)+x(prows+2)*b(1)+x(prows+3)*b(2);
  double sum=0;
  for(int i=0;i<prows;i++){
    double r=sqrt(pow(b(0)-P(i,0),2)+pow(b(1)-P(i,1),2)+pow(b(2)-P(i,2),2));
    sum=sum+x(i)*r;
  }
  return p+sum;
  
}
int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argc>1?argv[1]: "./bunny.off",V,F);
  int vrows=V.rows();
  V.resize(vrows,3);
  int frows=F.rows();
  F.resize(frows,3);
  double maxx=V(0,0); double minx=V(0,0);
  double maxy=V(0,1); double miny=V(0,1);
  double maxz=V(0,2); double minz=V(0,2);
  for(int i=1;i<vrows;i++){
    if (V(i,0)>maxx) maxx=V(i,0);
    if (V(i,0)<minx) minx=V(i,0);
    if (V(i,1)>maxy) maxy=V(i,1);
    if (V(i,1)<miny) miny=V(i,1);
    if (V(i,2)>maxz) maxz=V(i,2);
    if (V(i,2)<minz) minz=V(i,2);
   }
  cout<<"Showing BB"<<endl;
  cout<<maxx<<"  "<<maxy<<"  "<<maxz<<"  "<<endl;
  cout<<minx<<"  "<<miny<<"  "<<minz<<"  "<<endl;
  double bb=sqrt(pow(maxx-minx,2)+pow(maxy-miny,2)+pow(maxz-minz,2));
  cout<<bb<<endl;
  double r=0.01*bb;
  int rpoints=1000;
  Eigen::VectorXi I(rpoints);
  // cout<<"Showing Index"<<endl;
  // srand(0);
  // for(int i=0;i<rpoints;i++){
  //    I(i)=(rand()%vrows);
  //    //cout<<I(i)<<endl;
  //   }
  Eigen::MatrixXd R;
  R.resize(rpoints,3);
  // for(int i=0;i<rpoints;i++){
  //   R.row(i)=V.row(30*i);
  // }
    int prows=R.rows();  
    MatrixXd FN;
    FN.resize(vrows,3);
    cout<<"Calculating Normals"<<endl;
    igl::per_vertex_normals(V,F,FN);
    MatrixXd N;
    // for (int i=0;i<prows;i++){
    //     N.row(i)=FN.row(30*i);
    // }
    //int prows=vrows;
    MatrixXd A(2*prows,2*prows);
    double maximum;
    MatrixXd GV; 
    Eigen::RowVector3i res;
    Eigen::VectorXd S;
    do
    
    {
      cout<<"Number of vertices are: "<<rpoints<<endl;
      R.resize(rpoints,3);
      for(int i=0;i<rpoints;i++){
          R.row(i)=V.row(2*i);
        }
      int prows=R.rows(); 

      N.resize(prows,3);
      for (int i=0;i<prows;i++){
          N.row(i)=FN.row(2*i);
        }

      cout<<"Number of vertices are: "<<prows<<endl;
      MatrixXd P1(prows,3);
      P1.resize(prows,3);
      for (int i=0;i<prows-1;i=i+2){
        P1.row(i)=R.row(i)+r*N.row(i);
        P1.row(i+1)=R.row(i)-r*N.row(i);
      }
    MatrixXd P2(prows,3);
    P2.resize(prows,3);
    for (int i=1;i<prows;i=i+2){
      P2.row(i-1)=R.row(i)+r*N.row(i);
      P2.row(i)=R.row(i)-r*N.row(i);
     }
    MatrixXd P(2*prows,3);
    P.resize(2*prows,3);
    P<<P1,P2;
    cout<<P.rows()<<" "<<P.cols()<<endl;
   
    A.resize(2*prows,2*prows);
    
   cout<<"Calculating Euclidean Matrix"<<endl;

  for(int i=0;i<2*prows;i++){
    for(int j=i;j<2*prows;j++){
     //A(i,j)=sqrt(P.row(i).dot(P.row(i))+P.row(j).dot(P.row(j))-2*P.row(i).dot(P.row(j)));
      if(i==j) A(i,j)=0;
      A(i,j)=sqrt(pow(P(i,0)-P(j,0),2)+pow(P(i,1)-P(j,1),2)+pow(P(i,2)-P(j,2),2));
      A(j,i)=A(i,j);
    }
  }
  Eigen::MatrixXd P_(2*prows,4);
  for(int i=0;i<2*prows;i++){
    P_(i,0)=1.0;
    P_(i,1)=P(i,0);
    P_(i,2)=P(i,1);
    P_(i,3)=P(i,2);
  }
  Eigen::MatrixXd Pt;
  Pt=P_.transpose();
  Eigen::MatrixXd Zero=MatrixXd::Constant(4,4,0);
  Eigen::MatrixXd Fl(2*prows+4,2*prows+4);
  Fl<<A,P_,Pt,Zero;
  Eigen::MatrixXd b(2*prows+4,1);
  for(int i=0;i<2*prows-1;i=i+2){
    b(i,0)=0.05;
    b(i+1,0)=-0.05;
  }
  for(int i=2*prows;i<2*prows+4;i++){
    b(i,0)=0;
  }
  cout<<"Solving the system"<<endl;
  Eigen::VectorXd x=Fl.partialPivLu().solve(b);
  //for(int i=0;i<prows;i++){
    //cout<<x(i)<<endl;
    //}
  const int s = 50;
  // create grid
  
  cout<<"Computing Grid"<<endl;
  igl::voxel_grid(R,0,s,1,GV,res);
  cout<<GV.rows()<<"   "<<GV.cols()<<endl;
  cout<<res.rows()<<"   "<<res.cols()<<endl;

  S.resize(GV.rows());
  
  cout<<"Computing Interpolant values"<<endl;
  for(int i=0;i<GV.rows();i++){
    VectorXd u=GV.row(i);
    S(i)=interpolant(u,x,P);
    cout<<S(i)<<endl;
   }
  //  FILE *S1;
  //  S1=fopen("voxelgrid.txt","w");
  //  fwrite(&S,sizeof(double),1,S1);
  //  VectorXd p=P.row(2);
  //  cout<<"test "<<interpolant(p,x,P)<<endl;
   cout<<"Check"<<endl;
  std::vector<double> e;
  for(int i=0;i<vrows;i++){
    double residual=abs(0.05-interpolant(V.row(i),x,P));
    e.push_back(residual);
  }

  // for(int i=0; i<vrows;i++){
  //     cout<<e[i]<<endl;
  //     if(e[i]>0.23){
  //       prows=prows+1;
  //       R.resize(prows,3);
  //       R.row(prows-1)=V.row(i);
  //       N.resize(prows,3);
  //       N.row(prows-1)=FN.row(i);
  //     }
  // }
  maximum= *max_element(e.begin(),e.end()); 
  cout<<"Maximum of error is : "<<maximum<<endl;
  // if(prows>1500){
  //   break;
  // }
  rpoints=rpoints+200;
  if (2*rpoints>vrows){
    break;
  }
  
  

  }while(maximum>0.05);

  MatrixXd SV;
  MatrixXi SF;
  cout<<"Marching cubes"<<endl;
  igl::marching_cubes(S,GV,res(0),res(1),res(2),0,SV,SF);
  cout<<SV.rows()<<"  "<<SV.cols()<<endl;
  cout<<SF.rows()<<"  "<<SF.cols()<<endl;
  cout<<res(0)<<" "<<res(1)<<" "<<res(2)<<endl;
  // double dist;
  // igl::hausdorff(V,F,SV,SF,dist);
  // cout<<"Hausdorff distance is: "<<dist<<endl;
  //igl::writeOBJ("output.obj",SV,SF);
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(SV,SF);
  viewer.launch();
  
  
  
   
}
  
  
