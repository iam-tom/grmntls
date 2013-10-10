#include "garmin.h"
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <Eigen/Eigen>


bool DEBUG=false;

double interpolate_split(Eigen::VectorXd& dist,Eigen::VectorXd& time,double split_dist)
{
  double split_time=0.0;

  Eigen::VectorXd temp;
  temp=dist.array()-split_dist;
  temp=abs(temp.array());

  // get the elements for interpolation
  int min_index;
  double xa,xb,ya,yb;
  xa=temp.minCoeff(&min_index);
  if(xa >0)
  {
      xa=dist[min_index];
      ya=time[min_index];
      xb=dist[min_index+1];
      yb=time[min_index+1];
  }
   
  else if (xa <0)
  {
      xa=dist[min_index-1];
      ya=time[min_index-1];
      xb=dist[min_index];
      yb=time[min_index];
  }
  else
  {
    // no need for interpolation
    std::cerr<<"no int"<<std::endl;
    split_time=ya;
    return split_time;
  }


  //std::cout<<" xa:"<<xa<<" ya: "<<ya<<"\n xb:"<<xb<<" yb: "<<yb<<std::endl;
  split_time=ya+( ((yb-ya) / (xb-xa)) * (split_dist-xa));
  return split_time;

}
void calc_splits(Eigen::VectorXd& t,Eigen::VectorXd& d)
{


  //eliminate last entries

  std::cout<<" Aggregated Kilometer splits"<<std::endl;
  std::cout<<"[m]          "<<"[min] [sec]"<<std::endl;
  std::cout<<"---------------------------"<<std::endl;
  //reduce time
  t=t.array()-t(1,0);
  int index;
  d=d.head(d.size()-1);
  double max_dist=d[d.size()-1];

  double dist=1000;
  while(dist<max_dist)
  {

  double split_t2=interpolate_split(d,t,dist);
  std::cout<<dist<<" ......... "<<floor((split_t2)/60)<<":"<<floor((split_t2-floor(split_t2))*60)<<" "<<std::endl;
  //std::cout<<"-----------------------------------------------------------------"<<std::endl;
  
  dist+=1000;
  }
}
void extract_laps(garmin_data * gdata,std::vector<double>& t)
{

  // list processing
  if (gdata->type==1)
  {
    garmin_list_node * node;
    garmin_list * list;

    // convert data to list format
    gdata=garmin_list_data(gdata,1);


    if(gdata->type==data_Dlist)
    {
      list=static_cast<garmin_list* >(gdata->data);
    }
    else
    {
      std::cout<<"data is: "<<gdata->type<<std::endl;
    }

    // initialize vectors
    
    for ( node=list->head;node!=NULL;node=node->next)
    {
      //std::cout<<"Filetype: "<<node->data->type<<std::endl;
      switch (node->data->type)
      {
        case 1015:
          {
            D1015 * d1015;
            garmin_data * lap;
            lap=node->data;
            d1015=static_cast<D1015*>(lap->data);
            t.push_back(static_cast<double>(d1015->start_time));
            break;
          }
        default:
          {
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
void extract_points(garmin_data * gdata,std::vector<double>& d,std::vector<double>& t)
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
    else
    {
      std::cout<<"data is: "<<gdata->type<<std::endl;
    }

    // initialize vectors
    
    for ( node=list->head;node!=NULL;node=node->next)
    {
      //std::cout<<"Filetype: "<<node->data->type<<std::endl;
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
            //std::cout<<"processing for >>"<<node->data->type<<"<< not implemented"<<std::endl;
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
  std::vector<double>d_total;
  std::vector<double>t_total;
  std::vector<double>t_lap;
  extract_points(gdata,d_total,t_total);
  extract_laps(gdata,t_lap);


  Eigen::Map<Eigen::VectorXd> d_total_map(d_total.data(), d_total.size());
  Eigen::Map<Eigen::VectorXd> t_total_map(t_total.data(), t_total.size());
  Eigen::VectorXd d_total_vec(d_total.size());
  Eigen::VectorXd t_total_vec(d_total.size());
  d_total_vec=d_total_map;
  t_total_vec=t_total_map;
  calc_splits(t_total_vec,d_total_vec);
  //garmin_print_data(gdata,ofile,1);

  return 0;
}


