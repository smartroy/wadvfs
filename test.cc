#include "Core.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>

std::vector<std::ifstream* > files;
std::ofstream output("power_out_test");
std::ofstream output_slack("dynamic_slack");
std::ofstream output_tmpr("tmpr_c1");
std::ofstream output_rel("rel_c1");
std::ofstream output_ipc("ipc_c1");
std::ofstream output_dvfs("output_dvfs");
int main(int argc, char* argv[])
{
	Core test_core;
	//test_core.dis_temperature();
	std::vector<Task> tasks;
	std::ifstream workload;
	double sys_time=0;
	if(argc<3){
		std::cerr<<"Usage:"<<argv[0]<<" workload_file dvfs_mode"<<std::endl;
		return 1;
	}
	long lcm_total=1;
	
	workload.open(argv[1]);
	std::string workload_string;
	int period_ratio=std::atoi(argv[3]);

	while(getline(workload,workload_string)){
		Task* new_task=new Task;
		std::istringstream workload_line(workload_string);
		workload_line>>new_task->name>>new_task->index>>new_task->period>>new_task->wcet>>new_task->arrival>>new_task->actual>>new_task->mean>>new_task->high;
		new_task->period=new_task->period/10*period_ratio;
		new_task->deadline=new_task->period;
		new_task->name="power_file/"+new_task->name+".txt";
		new_task->exe=0;
		new_task->read_ratio=1;
		std::ifstream* new_file=new std::ifstream;
		new_file->open(new_task->name.c_str());
		files.push_back(new_file);
		tasks.push_back(*new_task);
		lcm_total=get_lcm(new_task->period,lcm_total);
	}
	//std::string test_string;
	//getline(*files[0],test_string);
	//std::cout<<test_string<<std::endl;
	for(int i=0;i<tasks.size();++i){
		std::cout<<tasks[i].name<<" "<<tasks[i].index<<" "<<tasks[i].period<<" "<<tasks[i].wcet<<" "<<tasks[i].deadline<<std::endl;
	}
	std::string dvfs_mode=argv[2];
	std::cout<<argv[1]<<std::endl;
	if(dvfs_mode=="pdvfs"){
		std::cout<<"uisng pdvfs"<<std::endl;
		test_core.set_dvfs_mode(true);
	}
	else{
		std::cout<<"using other dvfs"<<std::endl;
		test_core.set_dvfs_mode(false);
	}
	std::cout<<"period scale ratio "<<period_ratio<<std::endl;
	test_core.gen_lcm(tasks,1.5);
	//return 0;
	//test_core.run_lcm(lcm_total*3,sys_time);
	int cycle_count=0;
	output_rel<<argv[1]<<std::endl;
	while(test_core.rel_total>0.999999){
		test_core.run_lcm(lcm_total*3,sys_time);
		//std::cout<<"new round\n"<<std::endl;
		sys_time+=lcm_total*3*TIME_STEP;
		cycle_count+=1;
	}
	output<<"expired slack "<<test_core.expired_slack_total/cycle_count<<std::endl;
	//Task test_t={"basicmath.txt",1,10,5,0,0};
	//test_core.gen_lcm(&test_t,1);
	//test_core.run_lcm();
	//if (!test_t.file){
	//test_t.file.open(test_t.name.c_str());
	//}
	/*int data;
	std::string s;
	getline(test_t.file,s);
	std::istringstream line(s);
	std::string t_name;
	int data2;
	line>>data>>t_name>>data2;
	std::cout<<"test start "<<data<<" "<<t_name<<" "<<data2<<std::endl;
	int counter=4;
	while(getline(test_t.file,s)&&counter>0){
		//getline(test_t.file,s);
		std::istringstream line(s);
		line>>data>>t_name>>data2;
		std::cout<<"test while "<<data<<" "<<t_name<<" "<<data2<<std::endl;
		counter-=1;
	}
	while(getline(test_t.file,s)){
		
		std::istringstream line(s);
		line>>data>>t_name>>data2;
		std::cout<<"test end "<<data<<" "<<t_name<<" "<<data2<<std::endl;
		
	}*/
	//test_core.print_RC();
	/* code */
	return 0;

}