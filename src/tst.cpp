#include "garmin.h"
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <Eigen/Eigen>


bool DEBUG=false;

void calc_splits(Eigen::VectorXd& t,Eigen::VectorXd& d)
{

  //reduce time
  t=t.array()-t(1,0);

  Eigen::VectorXd temp;
  temp=d.array()-5000;
  temp=abs(temp.array());

  int min_index;
  double min;
  min=temp.minCoeff(&min_index);
  //std::cout<<"min= "<<min<<" at "<<min_index<<std::endl;

  std::cout<<"5000M split "<< t[min_index]/60<<std::endl;




}
void process_gdata(garmin_data * gdata,std::vector<double>& d,std::vector<double>& t)
{

  // list processing
  if (gdata->type==1)
  {
    garmin_list_node * node;
    garmin_list * list;

    // convert data to list format
    gdata=garmin_list_data(gdata,2);


    if(gdata->type==data_Dlist)
    {
      list=static_cast<garmin_list* >(gdata->data);
    }

    // initialize vectors
    
    for ( node=list->head;node!=NULL;node=node->next)
    {
      if (DEBUG)std::cout<<"Filetype: "<<node->data->type<<std::endl;
      switch (node->data->type)
      {
        case 304:
          {
            D304 * d304;
            garmin_data * point;
            point=node->data;
            d304=static_cast<D304*>(point->data);
            d.push_back((double)d304->distance);
            t.push_back(static_cast<double>(d304->time));
            break;
          }
        default:
          {
            std::cout<<"processing for >>"<<node->data->type<<"<< not implemented"<<std::endl;
            break;
          }
      }
    }


  }
  // when garmin data is no list
  else
  {
    std::cout<<"no list - exit"<<std::endl;
    exit(1);
  }

}

int main(int argc, const char *argv[])
{
 
  garmin_data* gdata;
  FILE* ofile;

  boost::filesystem::path input_path(argv[1]);
  boost::filesystem::path output_path;
  output_path=input_path;
  output_path.replace_extension(".dmp");
  std::cout<<"OUTPUT PATH: "<<output_path.c_str()<<std::endl;


  ofile=fopen(output_path.c_str(),"w");

  gdata=garmin_load(argv[1]);
  std::vector<double>d;
  std::vector<double>t;
  process_gdata(gdata,d,t);

  if (DEBUG)
  {
  for(int i =0;i<d.size();i++)
  {
    std::cout<<d[i]<<" " <<t[i]-t[0]<<std::endl;
  }
  }





  Eigen::Map<Eigen::VectorXd> d_map(d.data(), d.size());
  Eigen::Map<Eigen::VectorXd> t_map(t.data(), t.size());
  Eigen::VectorXd d_vec(d.size());
  Eigen::VectorXd t_vec(d.size());
  d_vec=d_map;
  t_vec=t_map;
  calc_splits(t_vec,d_vec);
  //garmin_print_data(gdata,ofile,1);

  return 0;
}


