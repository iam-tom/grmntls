#include "garmin.h"
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <Eigen/Eigen>


bool DEBUG=false;

double interpolate_split(Eigen::VectorXd& x,Eigen::VectorXd& y,double x0)
{
  double y0=0.0;

  Eigen::VectorXd temp;
  temp=x.array()-x0;
  temp=abs(temp.array());

  // get the elements for interpolation
  int min_index;
  double xa,xb,ya,yb;
  xa=temp.minCoeff(&min_index);
  if(xa >0)
  {
      xa=x[min_index];
      ya=y[min_index];
      xb=x[min_index+1];
      yb=y[min_index+1];
  }
   
  else if (xa <0)
  {
      xa=x[min_index-1];
      ya=y[min_index-1];
      xb=x[min_index];
      yb=y[min_index];
  }
  else
  {
    // no need for interpolation
    y0=y[min_index];
    return y0;
  }


  //std::cout<<" xa:"<<xa<<" ya: "<<ya<<"\n xb:"<<xb<<" yb: "<<yb<<std::endl;
  y0=ya+( ((yb-ya) / (xb-xa)) * (x0-xa));
  return y0;

}

void sec2hminsec(double& dec,std::vector<int>& h_min_sec)
{
  h_min_sec.resize(3);
  h_min_sec[0]=floor(dec/3600);
  h_min_sec[1]=floor((dec-h_min_sec[0]*3600)/60);
  h_min_sec[2]=floor(dec-(h_min_sec[0]*3600) - (h_min_sec[1]*60) );

  return;
}

void analyze_laps(Eigen::VectorXd t,Eigen::VectorXd d,Eigen::VectorXd l)
{

  Eigen::VectorXd lap_dist_vec;
  lap_dist_vec.resize(l.size());
  t=t.array()-t[0];
  l=l.array()-l[0];

  for(int no_lap=0;no_lap<l.size();no_lap++)
  {
  double lap_dist=interpolate_split(t,d,l[no_lap]);
  lap_dist_vec[no_lap]=lap_dist;

  if (no_lap>0)
  {
    double delta_t=l[no_lap]-l[no_lap-1];
    std::vector<int>hex;
    sec2hminsec(delta_t,hex);
    double delta_d=lap_dist_vec[no_lap]-lap_dist_vec[no_lap-1];


    std::cout<<"time for lap "<<no_lap<<" = "<<hex[1]<<","<<hex[2]<<" distance = "<<delta_d<<std::endl;
  }
  }

}

void fastest_xm(Eigen::VectorXd t, Eigen::VectorXd d,int x)
{
  // fastest km between   >------|.....|------->
  //                      d0    da     db    dmax
  double dmax=d[d.size()-2];
  std::cout<<"DMAX"<<dmax<<std::endl;

  double da,db,da_min,db_min;//=d[0]+1000;
  int i=0;
  double dt_min=std::numeric_limits<double>::max();
  da=0;
  db=0;
  double ta,tb,dt;

  for(int i=0;i<d.size();++i)
  {
      da=d[i];
      db=d[i]+x;
      if(db<0 || db>1e13) continue;
      if(db>dmax)  break;
      ta=t[i];
      tb=interpolate_split(d,t,db);
      dt=tb-ta;

      if(dt<dt_min)
      {
        dt_min=dt;
        da_min=da;
        db_min=db;
      }
  }

    std::vector<int> hex;
    sec2hminsec(dt_min,hex);
    std::cout<<"Fastest "<<x<<" m"<<std::endl;
    std::cout<<da_min<<" --> "<<db_min<<" | "<<hex[1]<<"min "<<hex[2]<<"s"<<std::endl; 



}
void fastest_km(Eigen::VectorXd t, Eigen::VectorXd d)
{
  // fastest km between   >------|.....|------->
  //                      d0    da     db    dmax
  double dmax=d[d.size()-2];
  std::cout<<"DMAX"<<dmax<<std::endl;

  double da,db,da_min,db_min;//=d[0]+1000;
  int i=0;
  double dt_min=std::numeric_limits<double>::max();
  da=0;
  db=0;
  double ta,tb,dt;

  for(int i=0;i<d.size();++i)
  {
      da=d[i];
      db=d[i]+1000.0;
      if(db<0 || db>1e13) continue;
      if(db>dmax)  break;
      ta=t[i];
      tb=interpolate_split(d,t,db);
      dt=tb-ta;

      if(dt<dt_min)
      {
        dt_min=dt;
        da_min=da;
        db_min=db;
      }
  }

    std::vector<int> hex;
    sec2hminsec(dt_min,hex);
    std::cout<<"Fastest km"<<std::endl;
    std::cout<<da_min<<" --> "<<db_min<<" | "<<hex[1]<<"min "<<hex[2]<<"s"<<std::endl; 



}
  

void calc_splits(Eigen::VectorXd t,Eigen::VectorXd d)
{


  //eliminate last entries

  std::cout<<" Aggregated Kilometer splits"<<std::endl;
  std::cout<<"[m]          "<<"[h] [min] [sec]"<<std::endl;
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
  //std::cout<<dist<<" ......... "<<floor((split_t2)/60)<<":"<<floor((split_t2-floor(split_t2))*60)<<" "<<std::endl;
  std::vector<int>time_hex;
  sec2hminsec(split_t2,time_hex);
  std::cout<<dist<<" ......... "<<time_hex [0]<<" "<<time_hex [1]<<":"<<time_hex[2]<<" "<<std::endl;




  //std::cout<<"-----------------------------------------------------------------"<<std::endl;
  
  dist+=1000;
  }
}
void extract_runs(garmin_data * gdata,std::vector<double>& t)
{

  // list processing
  if (gdata->type==1)
  {
    garmin_list_node * node;
    garmin_list * list;

    // convert data to list format
    gdata=garmin_list_data(gdata,0);


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
        case 1009:
          {
            D1009 * d1009;
            garmin_data * run;
            run=node->data;
            d1009=static_cast<D1009*>(run->data);
            break;
          }
        default:
          {
            break;
          }
      }
    }


  }
  // when garmin data is no list but single element
  else
  {
      switch (gdata->type)
      {
        case 1009:
          {
            D1009 * d1009;
            d1009=static_cast<D1009*>(gdata->data);


          }
        default:
          {
            break;
          }
      }


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
      switch (gdata->type)
      {
        case 1015:
          {
            D1015 * d1015;
            garmin_data * lap;
            d1015=static_cast<D1015*>(gdata->data);
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
void extract_points(garmin_data * gdata,std::vector<double>& d,std::vector<double>& t)
{
  int i=1;
  i+=1;

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
      std::cout<<i<<"data is: "<<gdata->type<<std::endl;
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


  //ofile=fopen(output_path.c_str(),"w");

  gdata=garmin_load(argv[1]);
  std::vector<double>d_total;
  std::vector<double>t_total;
  std::vector<double>t_lap;
  extract_points(gdata,d_total,t_total);
  extract_laps(gdata,t_lap);


  Eigen::Map<Eigen::VectorXd> d_total_map(d_total.data(), d_total.size());
  Eigen::Map<Eigen::VectorXd> t_total_map(t_total.data(), t_total.size());
  Eigen::Map<Eigen::VectorXd> t_lap_map(t_lap.data(), t_lap.size());

  Eigen::VectorXd d_total_vec(d_total.size());
  Eigen::VectorXd t_total_vec(d_total.size());
  Eigen::VectorXd t_lap_vec(t_lap.size());

  d_total_vec=d_total_map;
  t_total_vec=t_total_map;
  t_lap_vec=t_lap_map;

  calc_splits(t_total_vec,d_total_vec);
  analyze_laps(t_total_vec,d_total_vec,t_lap_vec);
  fastest_km(t_total_vec,d_total_vec);
  fastest_xm(t_total_vec,d_total_vec,5000);
  fastest_xm(t_total_vec,d_total_vec,10000);
  fastest_xm(t_total_vec,d_total_vec,500);

  //garmin_print_data(gdata,ofile,1);

  return 0;
}


