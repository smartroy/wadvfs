//when task expire, static slack from this task renew. in handle_expire, after deal with expire_dynamic, reset the static of this
//task by add a new variable hold all static, decrease if no dyanmic slack available in use_slack, when a task renew, 
//s_h+=(1-all_static_h/s_static_h)*static_of_task*(s_static_h/total_static) s_l+=...

//use slack that about to expire, dynamic: when expire-slack/step_use>=core_time, add all slack to s_left. static: create a 
// deadline_array[],store the deadline of current invoke, when unused_static_l(h)/static_total*static_of_task_i will not be used by
//expire, move all to s_l

#include "Core.h"
#include <iostream>
#include <deque>
#include <cmath>
#include <iomanip>
extern "C"{
	#include "flp.h"
	#include "temperature.h"
	#include "tilts.h"
}
#include <algorithm>
#include <chrono>
#include <random>
//#include "MersenneTwister.o"
//#include <time.h>
const double SCALE_EM(0.3);
const double SCALE_BD(0.18);
const double K_BOL(8.62e-5);
const double EA_EM(0.9);
const double EXP_EM(2.0);
const double EXP_BD(2.0);
Core::Core(){
	set_temperature();
	gen_flp();
}
extern std::vector<std::ifstream* > files;
extern std::ofstream output;
extern std::ofstream output_slack;
extern std::ofstream output_tmpr;
extern std::ofstream output_rel;
extern std::ofstream output_ipc;
extern std::ofstream output_dvfs;
void Core::set_temperature(){
	for(int i=0;i<blocks;i++){
		temperature[i]=T_INIT+273.15;
		power[i]=0.0;
		power_cal[i]=0.0;
		tmpr[i]=T_INIT+273.15;
		area_ratio[i]=1;
	}
	for(int i=blocks;i<EXTRA+blocks;i++){
		power_cal[i]=0.0;
	}
	for(int i=blocks;i<NL*blocks+EXTRA;i++){
		tmpr[i]=T_INIT+273.15;
	}
	
	spd[LV]=6;
	spd[MV]=6;
	spd[HV]=10;
	p_ratio[LV]=0.41;//0.41;//0.32;0.54;0.68;
	p_ratio[MV]=0.41;
	p_ratio[HV]=1.0;
	rel_em=1;
	rel_bd=1;
	rel_cycle=1;
	rel_total=1;
	expired_slack_total=0;

}

void Core::set_temperature(float *temp_array){
	for(int i=0;i<blocks;++i){
		temperature[i]=temp_array[i];
	}
}

void Core::dis_temperature(){
	for(int i=0;i<blocks;i++){
		std::cout<<temperature[i]<<std::endl;
	}
}

void Core::gen_flp(){
	flp=read_flp("90nm.flp",FALSE);
	config=default_thermal_config();
	rc_model=alloc_RC_model(&config,flp);
	populate_R_model(rc_model,flp);
	populate_C_model(rc_model,flp);
	initialize_tilts(rc_model,flp,rc_model->config->sampling_intvl,Aj,Bx);

}
void Core::print_RC(){
	debug_print_model(rc_model);
}
void Core::gen_lcm(std::vector<Task> t_list,float threshold){
	
	//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	//std::default_random_engine generator(seed);
  	
	ipc_threshold=threshold;
	core_utilization=0;
	s_static=0;
	static_l_aet=0;
	static_h_aet=0;
	static_l_wcet=0;
	static_h_wcet=0;
	static_left=0;
	s_h_aet=0;
	s_l_aet=0;
	s_h_wcet=0;
	s_l_wcet=0;

	s_static_h=0;
	s_static_l=0;
	exe_h=0;
	exe_l=0;
	core_lcm=t_list[0].period;
	int reserve_voltage=LV;
	//init_genrand(time(NULL))
	task_count=t_list.size();
	for(int i=0;i<t_list.size();++i){
		std::cout<<"lcm of "<<core_lcm<<" "<<t_list[i].period<<std::endl;
		core_lcm=get_lcm(core_lcm,t_list[i].period);
		output<<"core lcm "<<core_lcm<<std::endl;
		core_utilization+=float(t_list[i].wcet)/float(t_list[i].period*10);
		//core_high+=float(t_list[i].wcet)*t_list[i].high;
		//core_wcet+=t_list[i].wcet;
		task_utilization[t_list[i].index]=float(t_list[i].wcet)/float(t_list[i].period*10);
		Defer_task* temp_defer=new Defer_task;
		temp_defer->index=t_list[i].index;
		temp_defer->left=t_list[i].wcet;
		temp_defer->deadline=t_list[i].deadline;
		temp_defer->utilization=float(t_list[i].wcet)/float(t_list[i].period*10);
		defer_list.push_back(*temp_defer);
	}
	wc_utilization=core_utilization;
	task_list.clear();
	for(int i=0;i<t_list.size();++i){
		
		for(int j=0;j<core_lcm/t_list[i].period;j++){
			/*std::normal_distribution<double> distribution (t_list[i].mean,0.2);
			double exe_rand=0;
			while(exe_rand<=0|| exe_rand>1.0){
				exe_rand=distribution(generator);
			}*/
			Task *temp_task=copy_task(t_list[i],j,t_list[i].actual);
			task_list.push_back(*temp_task);
		}
		s_static+=float(t_list[i].wcet)*(1/core_utilization-1);
		static_array[t_list[i].index]=float(t_list[i].wcet)*(1/core_utilization-1);
		s_h_aet+=float(t_list[i].wcet)*t_list[i].mean*t_list[i].high/spd[reserve_voltage]*(spd[HV]-spd[reserve_voltage]);
		s_l_aet+=float(t_list[i].wcet)*t_list[i].mean*(1.0 - t_list[i].high)/spd[reserve_voltage]*(spd[HV]-spd[reserve_voltage]);
		s_h_wcet+=float(t_list[i].wcet)*(1.0 - t_list[i].mean)*t_list[i].high/spd[reserve_voltage]*(spd[HV]-spd[reserve_voltage]);
		s_l_wcet+=float(t_list[i].wcet)*(1.0 - t_list[i].mean)*(1 - t_list[i].high)/spd[reserve_voltage]*(spd[HV]-spd[reserve_voltage]);
		exe_h+=float(t_list[i].wcet)*t_list[i].high;
		exe_l+=float(t_list[i].wcet)*(1.0 - t_list[i].high);
		dynamic_array[t_list[i].index]=0;
		output<<" "<<s_h_aet;

		output_slack<<static_array[t_list[i].index]<<" ";
	}
	static_total=s_static;
	output_slack<<std::endl;
	output<<std::endl;
	for(int i=0;i<task_list.size();++i){
		output<<task_list[i].name<<" "<<task_list[i].arrival<<" "<<task_list[i].actual<<std::endl;
	}
	output<<"utilization "<<core_utilization<<" slack assign: total "<<s_static;
	if(dvfs_tradition){
		cal_static();
	}
	else{
		cal_static();
	}
	dynamic_reset();
	std::sort(task_list.begin(),task_list.end(),cmp_arrival);
	
	output<<" high aet/wcet "<<s_h_aet<<"/"<<s_h_wcet<<" low aet/wcet "<<s_l_aet<<"/"<<s_l_wcet<<" high "<<s_h<<" low "<<s_l<<std::endl;
	
	//reset slack
	

}
void Core::run_lcm(long sys_wide_lcm,double sys_time){
	/*for(int i=0;i<task_list.size();++i){
		std::cout<<task_list[i].name<<" "<<task_list[i].period<<std::endl;
	}*/
	int task_pt=0;
	sys_lcm=sys_wide_lcm;
	core_time=0;
	int output_count=0;
	//std::string power_string;
	core_speed=spd[HV];
	int dvfs_count=0;
	unsigned seed =5;//std::chrono::system_clock::now().time_since_epoch().count();
  	std::default_random_engine generator(seed);
	IPC=0;
	IPC_dvfs=0;
	voltage=HV;
	std::cout<<"new round "<<task_list[0].arrival<<" "<<std::endl;
	int new_release=0;
	float out_ipc=0;
	while(core_time<sys_lcm){
		if(core_time%core_lcm==0){
			task_pt=0;
			if(!dvfs_tradition)
				dynamic_reset();
						
		}
		//output<<core_time<<" ";
		
		while(task_pt<task_list.size()&&core_time>=task_list[task_pt].arrival){
			//std::cout<<"adding task to exe "<<core_time<<task_list[0].actual<<std::endl;
			std::normal_distribution<double> distribution (task_list[task_pt].mean,0.2);
			double exe_rand=0;
			while(exe_rand<1e-5|| exe_rand>1.0){
				exe_rand=distribution(generator);
			}
				
			Task* temp_task=copy_task(task_list[task_pt],0,int(float(task_list[task_pt].wcet)*exe_rand));
			temp_task->read_ratio=exe_rand;
			task_list[task_pt].arrival+=core_lcm;
			task_list[task_pt].deadline+=core_lcm;
			defer_array[temp_task->index].left=temp_task->wcet;
			defer_array[temp_task->index].deadline=temp_task->deadline;
			exe_list.push_back(*temp_task);
			core_utilization-=task_utilization[temp_task->index];
			task_utilization[temp_task->index]=float(temp_task->wcet)/float(temp_task->period*10);
			core_utilization+=float(temp_task->wcet)/float(temp_task->period*10);
			
			//task_list.pop_front();
			std::sort(exe_list.begin(),exe_list.end(),cmp_deadline);
			task_pt+=1;
			new_release=1;
		}
		if(dvfs_tradition&&new_release==1){
			//pdvfs_defer();
			set_speed();
			new_release=0;
		}
		if(dvfs_count>=dvfs_step&&(!dvfs_tradition)){
			IPC_dvfs/=dvfs_count;
			dvfs_count=0;
			handle_expire();
			if(new_dynamic>0){
				output_slack<<"new dynamic "<<new_dynamic<<" "<<dynamic_h_aet+static_h_aet<<"/"<<dynamic_h_wcet+static_h_wcet<<"/"<<dynamic_l_aet+static_l_aet<<"/"<<dynamic_l_wcet+static_l_wcet<<"/"<<dynamic_left+static_left<<std::endl;
				alloc_dynamic();
				output_slack<<"new dynamic2 "<<new_dynamic<<" "<<dynamic_h_aet+static_h_aet<<"/"<<dynamic_h_wcet+static_h_wcet<<"/"<<dynamic_l_aet+static_l_aet<<"/"<<dynamic_l_wcet+static_l_wcet<<"/"<<dynamic_left+static_left<<std::endl;
			}
			handle_about_expire();
			if(exe_list.size()>0){
				if(s_left>(spd[HV]-spd[LV])*dvfs_step){
					voltage=LV;
					core_speed=spd[voltage];
					s_left-=(spd[HV]-spd[LV])*dvfs_step;
					std::sort(slack_list.begin(),slack_list.end(),cmp_slack);				
					float s_used=(spd[HV]-spd[LV])*dvfs_step;
					use_slack_left(s_used);
					
				}
				else{
					if(IPC_dvfs>=ipc_threshold&&s_h>(spd[HV]-spd[LV])*dvfs_step){
						voltage=LV;
						core_speed=spd[voltage];
						s_h-=(spd[HV]-spd[LV])*dvfs_step;
						//use_slack function
						std::sort(slack_list.begin(),slack_list.end(),cmp_slack);				
						float s_used=(spd[HV]-spd[LV])*dvfs_step;
						use_slack_high(s_used);
						
					}
					else if(IPC_dvfs<ipc_threshold &&s_l>(spd[HV]-spd[LV])*dvfs_step){
						voltage=LV;
						core_speed=spd[voltage];
						s_l-=(spd[HV]-spd[LV])*dvfs_step;
						std::sort(slack_list.begin(),slack_list.end(),cmp_slack);				
						float s_used=(spd[HV]-spd[LV])*dvfs_step;
						use_slack_low(s_used);
						
					}
					else{
						voltage=HV;
						core_speed=spd[voltage];
					}
				}
			}
			else{
				voltage=LV;
				core_speed=spd[voltage];
			}
			out_ipc = IPC_dvfs;
			IPC_dvfs=0;
		}

		/*output_slack<<core_time<<" "<<"d_array ";
		for(int i=0;i<4;++i){
			output_slack<<dynamic_array[i]<<" ";
		}	
		output_slack<<" | "<<dynamic_h_aet_s<<"/"<<dynamic_h_wcet_s<<"/"<<dynamic_l_aet_s<<"/"<<dynamic_l_wcet_s<<"/"<<dynamic_left_s<<" | "<<static_h_total<<"/"<<static_l_total<<"/"<<static_left_total<<std::endl;*/
		core_time+=1;
		//output<<" IPC "<<IPC<<" speed "<<core_speed<<" slack h/l "<<s_h<<"/"<<s_l<<std::endl;
		
		if(exe_list.size()>0){
			//output<<exe_list.size()<<" ";
			IPC=0;		
			core_read(core_speed);
			
			//output<<" "<<exe_list.size()<<" ";

			IPC/=core_speed;
			IPC_dvfs+=IPC;
			//std::cout<<power_cal[3]<<" finishe reading "<<tmpr[3]<<std::endl;
			for(int i=0;i<blocks;++i){
				//std::cout<<" finishe reading "<<i<<core_speed<<p_ratio[voltage]<<std::endl;
				power_cal[i]=power_cal[i]/core_speed*p_ratio[voltage];
				
			}
		}
		output_dvfs<<IPC<<" "<<core_speed<<std::endl;
		output<<core_time<<" ";
		for(int i=0;i<task_count;++i){
			output<<task_utilization[i]<<" ";
			//std::cout<<power_cal[3]<<" finishe reading "<<tmpr[3]<<std::endl;
		}
		output<<core_utilization<<" "<<core_speed<<" "<<s_l<<" "<<s_h<<" "<<s_left<<std::endl;
		/*for(int i=0;i<defer_list.size();++i){
			output<<defer_array[i].left<<" "<<defer_array[i].deadline<<" ";
		}
		output<<std::endl;*/
		compute_temp_tilts(tmpr,power_cal,TIME_SLICE,N_BLOCKS,N_BLOCKS*NL+EXTRA,Aj,Bx);
		output_ipc<<IPC<<std::endl;
		cal_rel(sys_time+core_time*TIME_STEP);
		total_power=0;
		for(int i=0;i<blocks;++i){
			output_tmpr<<tmpr[i]<<" ";
			total_power+=power_cal[i];
			power_cal[i]=0;

		}
		output_tmpr<<std::endl;
		if(rel_total<0.999999){
			output_rel<<wc_utilization<<std::endl;
			output_rel<<sys_time+core_time*TIME_STEP<<" "<<total_power<<" "<<std::setprecision(9)<<std::scientific<<rel_total<<std::endl;
			break;
		}
		
		//output_slack<<core_time;
		if(slack_list.size()>0){
			//std::cout<<"printing dynamic"<<std::endl;
			std::deque<Slack>::iterator it_slack=slack_list.begin();
			while(it_slack!=slack_list.end()){
				//output_slack<<" "<<it_slack->expire<<" "<<it_slack->s;
				++it_slack;
			}
		}
		//output_slack<<std::endl;
		test_power=0;
		
		dvfs_count+=1;
		//apply dvfs

		

		
	}
	//reset arrival
	for(task_pt=0;task_pt<task_list.size();++task_pt){
		task_list[task_pt].arrival-=sys_lcm;
		task_list[task_pt].deadline-=sys_lcm;
	}
	
}
void Core::set_speed(){
	/*if(core_utilization>float(spd[MV])/float(spd[HV])){
		voltage=HV;
		core_speed=spd[voltage];
	}
	else */
	if(core_utilization>float(spd[LV])/float(spd[HV])){
		voltage=HV;
		core_speed=spd[voltage];
	}
	else{
		voltage=LV;
		core_speed=spd[voltage];
	}
	//	voltage=HV;
	//	core_speed=spd[voltage];
}
void Core::defer_set_speed(float ratio){
	if(ratio<float(spd[LV])/float(spd[HV])){
		voltage=LV;
		core_speed=spd[voltage];
	}
	else{
		voltage=HV;
		core_speed=spd[voltage];
	}
}
void Core::cal_rel(double rel_time){
	double alpha_em=0;
	double alpha_bd=0;
	double alpha_cycle=0;
	double bd_pre=1;
	double bd=1;
	double em_pre=1;
	double em=1;
	for(int i=0;i<blocks;++i){
		alpha_em=SCALE_EM/std::pow(1e6,1.1)*std::exp(EA_EM/(tmpr[i]*K_BOL));
		alpha_bd=SCALE_BD*16e7*std::pow(1/1.1,78+0.081*tmpr[i])*std::exp((0.759-66.8/tmpr[i]-8.37e-4*tmpr[i])/(K_BOL*tmpr[i]));
		bd_pre*=std::exp(-std::pow((rel_time - TIME_STEP)/alpha_bd,EXP_BD)*area_ratio[i]);
		bd*=std::exp(-std::pow(rel_time/alpha_bd,EXP_BD)*area_ratio[i]);
		em_pre*=std::exp(-std::pow((rel_time - TIME_STEP)/alpha_em,EXP_EM)*area_ratio[i]);
		em*=std::exp(-std::pow(rel_time/alpha_em,EXP_EM)*area_ratio[i]);

	}
	rel_total*=(1-(em_pre*bd_pre - em*bd));
}


void Core::core_read(int read_speed){
	std::string power_string;
	int count=0;
	int temp_power;
	float ipc_temp=0;
	if(exe_list[0].exe==0){
		files[exe_list[0].index]->close();
		files[exe_list[0].index]->open(exe_list[0].name.c_str());
	}
	int inner_speed=0;
	int inner_count=0;
	//std::cout<<"start reading "<<exe_list.size()<<std::endl;
	if (read_speed<=exe_list[0].actual){
		inner_speed=read_speed;
		inner_count=int(inner_speed/exe_list[0].read_ratio);
		exe_list[0].actual-=read_speed;
		exe_list[0].exe+=read_speed;
		read_speed=0;
	}
	else{
		inner_speed=exe_list[0].actual;
		inner_count=int(inner_speed/exe_list[0].read_ratio);
		read_speed-=inner_speed;
		exe_list[0].actual=0;
		exe_list[0].exe+=inner_speed;

	}
	while(count<inner_count){
		//std::cout<<exe_list[0].name<<std::endl;
		getline(*files[exe_list[0].index],power_string);
		std::istringstream line(power_string);
		//line>>temp_power;
		//test_power=temp_power;
		//output<<test_power<<" ";
		line>>ipc_temp>>power[0]>>power[1]>>power[2]>>power[3]>>power[4]>>power[5]>>power[6]>>power[7]>>power[8]>>power[9]>>power[10]>>power[11]>>power[12]>>power[13]>>power[14];
		for(int i=0;i<blocks;++i){
			power_cal[i]+=power[i]*inner_speed/inner_count;
		}
		IPC+=ipc_temp*inner_speed/inner_count;
		//exe_list[0].exe+=1;
		//exe_list[0].actual-=1;
		//defer_array[exe_list[0].index].left-=1;
		count+=1;
	}
	if(exe_list[0].actual==0){
		if(exe_list[0].deadline>=core_time + dvfs_step){
			//std::cout<<"new dynamic"<<std::endl;
			Slack* new_slack=new Slack;
			new_slack->expire=exe_list[0].deadline;
			new_slack->s=exe_list[0].wcet-exe_list[0].exe;
			new_slack->index=exe_list[0].index;
			slack_list.push_back(*new_slack);
			dynamic_array[exe_list[0].index]=exe_list[0].wcet-exe_list[0].exe;
			new_dynamic+=exe_list[0].wcet-exe_list[0].exe;
			
		}
		
		output<<"task finished "<<exe_list[0].exe<<" "<<exe_list[0].wcet<<" "<<exe_list[0].period;
		core_utilization-=task_utilization[exe_list[0].index];
		task_utilization[exe_list[0].index]=float(exe_list[0].exe)/float(exe_list[0].period*10);
		core_utilization+=task_utilization[exe_list[0].index];
		output<<" "<<task_utilization[exe_list[0].index]<<std::endl;
		//defer_array[exe_list[0].index].left=0;
		if(dvfs_tradition){
			//pdvfs_defer();
			set_speed();
		}
		exe_list.pop_front();
	}
	if(read_speed>0&&exe_list.size()>0){
		//exe_list.pop_front();
		core_read(read_speed);
	}
}

void Core::dynamic_reset(){
	dynamic_l_aet=0;
	dynamic_h_aet=0;
	dynamic_l_wcet=0;
	dynamic_h_wcet=0;
	dynamic_left=0;
	new_dynamic=0;
	dynamic_l_aet_s=0;
	dynamic_h_aet_s=0;
	dynamic_l_wcet_s=0;
	dynamic_h_wcet_s=0;
	dynamic_left_s=0;
	//expire_dynamic=0;
	
	s_l=static_l_aet + static_l_wcet;
	s_h=static_h_aet + static_h_wcet;
	s_left=static_left;

	static_l_total=static_l_aet + static_l_wcet;
	static_h_total=static_h_aet + static_h_wcet;
	static_left_total=static_left;
	
	slack_list.clear();
}
void Core::cal_static_tradition(){
	s_static_h=s_static*exe_h/(exe_h+exe_l);
	s_static_l=s_static*exe_l/(exe_h+exe_l);
}
void Core::cal_static(){
	if (s_h_aet>s_static){
		static_h_aet=s_static;
		s_static=0;
	}
	else{
		static_h_aet=s_h_aet;
		s_static-=s_h_aet;
		if(s_l_aet>s_static){
			static_l_aet=s_static;
			s_static=0;
		}
		else{
			static_l_aet=s_l_aet;
			s_static-=s_l_aet;
			if(s_h_wcet>s_static){
				static_h_wcet=s_static;
				s_static=0;
			}
			else{
				static_h_wcet=s_h_wcet;
				s_static-=s_h_wcet;
				if(s_l_wcet>s_static){
					static_l_wcet=s_static;
					s_static=0;
				}
				else{
					static_l_wcet=s_l_wcet;
					s_static-=s_l_wcet;
					static_left=s_static;
					s_static=0;
				}
			}
		}
	}
}
void Core::alloc_dynamic_tradition(){
	s_dynamic_h+=new_dynamic*exe_h/(exe_l+exe_h);
	s_dynamic_l+=new_dynamic*exe_h/(exe_h+exe_l);
	s_h+=new_dynamic*exe_h/(exe_l+exe_h);
	s_l+=new_dynamic*exe_l/(exe_l+exe_h);
	new_dynamic=0;
}
void Core::alloc_dynamic(){
	if(new_dynamic<s_h_aet-static_h_aet - dynamic_h_aet){
		dynamic_h_aet_s+=new_dynamic;
		dynamic_h_aet+=new_dynamic;
		s_h+=new_dynamic;
		new_dynamic=0;
	}
	else{
		if(s_h_aet>static_h_aet + dynamic_h_aet){
			s_h+=(s_h_aet - static_h_aet - dynamic_h_aet);
			dynamic_h_aet_s+=(s_h_aet - static_h_aet - dynamic_h_aet);
			new_dynamic-=(s_h_aet - static_h_aet - dynamic_h_aet);
			dynamic_h_aet+=(s_h_aet - static_h_aet - dynamic_h_aet);
			
			
		}
		if(new_dynamic<s_l_aet - static_l_aet - dynamic_l_aet){
			dynamic_l_aet_s+=new_dynamic;
			dynamic_l_aet+=new_dynamic;
			s_l+=new_dynamic;
			new_dynamic=0;
		}
		else{
			if(s_l_aet>static_l_aet + dynamic_l_aet){
				s_l+=(s_l_aet - static_l_aet - dynamic_l_aet);
				dynamic_l_aet_s+=(s_l_aet - static_l_aet - dynamic_l_aet);
				new_dynamic-=(s_l_aet - static_l_aet - dynamic_l_aet);
				dynamic_l_aet+=(s_l_aet - static_l_aet - dynamic_l_aet);
				
				
			}
			if(new_dynamic<s_h_wcet - static_h_wcet - dynamic_h_wcet){
				s_h+=new_dynamic;
				dynamic_h_wcet_s+=new_dynamic;
				dynamic_h_wcet+=new_dynamic;
				new_dynamic=0;
			}
			else{
				if(s_h_wcet > static_h_wcet + dynamic_h_wcet){
					s_h+= (s_h_wcet - static_h_wcet - dynamic_h_wcet);
					dynamic_h_wcet_s+=(s_h_wcet - static_h_wcet - dynamic_h_wcet);
					new_dynamic-=(s_h_wcet - static_h_wcet - dynamic_h_wcet);
					dynamic_h_wcet+=(s_h_wcet - static_h_wcet - dynamic_h_wcet);
					
					
				}	
				if(new_dynamic<s_l_wcet - static_l_wcet - dynamic_l_wcet){
					s_l+=new_dynamic;
					dynamic_l_wcet_s+=new_dynamic;
					dynamic_l_wcet+=new_dynamic;
					
					new_dynamic=0;
				}
				else{
					if(s_l_wcet>static_l_wcet + dynamic_l_wcet){
						s_l+=(s_l_wcet - static_l_wcet - dynamic_l_wcet);
						dynamic_l_wcet_s+=(s_l_wcet - static_l_wcet - dynamic_l_wcet);
						new_dynamic-=(s_l_wcet - static_l_wcet - dynamic_l_wcet);
						dynamic_l_wcet+=(s_l_wcet - static_l_wcet - dynamic_l_wcet);
						
						
					}
					dynamic_left+=new_dynamic;
					dynamic_left_s+=new_dynamic;
					s_left+=new_dynamic;
					new_dynamic=0;
				}
				
			}
		}
	}
}

void Core::pdvfs_defer(){
	std::deque<Defer_task>::iterator it=defer_list.begin();
	while(it!=defer_list.end()){
		it->left=defer_array[it->index].left;
		it->deadline=defer_array[it->index].deadline;
		++it;
	}
	std::sort(defer_list.begin(),defer_list.end(),defer_cmp_deadline_reverse);
	float U=wc_utilization;
	float x=0;
	float s=0;
	int size=defer_list.size();
	//output_slack<<core_time<<" ";
	for(int i=0;i<size;++i){
		if(defer_list[i].deadline>defer_list[size-1].deadline){
			U-=defer_list[i].utilization;
			x=(0>defer_list[i].left/10-(1-U)*float(defer_list[i].deadline - defer_list[size-1].deadline))?0:defer_list[i].left/10-(1-U)*float(defer_list[i].deadline - defer_list[size-1].deadline);
			U+=(defer_list[i].left/10 - x)/(defer_list[i].deadline - defer_list[size-1].deadline);
			s+=x;
		}
		else{
			x=defer_list[i].left/10;
			s+=x;
		}
		//output_slack<<"|"<<x<<" "<<defer_list[i].left/10<<" "<<float(defer_list[i].deadline)<<" |";
	}
	//output_slack<<s<<" "<<(s/float(defer_list[size-1].deadline-core_time))<<std::endl;
	if(defer_list[size-1].deadline>core_time){
		defer_set_speed(s/float(defer_list[size-1].deadline-core_time));
	}
	else{
		defer_set_speed(1.0);
	}
}
void Core::handle_expire(){
	float expire_dynamic=0;
	float expire_dynamic_total=0;
	if(slack_list.size()>0){
		std::deque<Slack>::iterator it=slack_list.begin();
		output_slack<<core_time<<" "<<s_h<<"/"<<s_l<<"/"<<s_left;
		while(it!=slack_list.end()){
			output_slack<<" | "<<it->expire<<" "<<it->s<<" "<<it->index;
			if(it->expire<core_time + dvfs_step ){
				output_slack<<" "<<" erase slack ";
				expired_slack_total+=(it->s/(spd[HV]-spd[LV]));
				expire_dynamic+=it->s;
				expire_dynamic_total+=dynamic_array[it->index];
				s_l+=(static_l_aet+static_l_wcet- static_l_total)*static_array[it->index]/static_total;
				s_h+=(static_h_aet+static_h_wcet- static_h_total)*static_array[it->index]/static_total;
				s_left+=(static_left- static_left_total)*static_array[it->index]/static_total;
				static_l_total+=(static_l_aet+static_l_wcet- static_l_total)*static_array[it->index]/static_total;
				static_h_total+=(static_h_aet+static_h_wcet- static_h_total)*static_array[it->index]/static_total;
				static_left_total+=(static_left- static_left_total)*static_array[it->index]/static_total;

				it=slack_list.erase(it);
			}
			else
				++it;
		}
		output_slack<<" | "<<expire_dynamic<<" "<<s_h<<"/"<<s_l<<"/"<<s_left<<std::endl;
		//output_slack<<expire_dynamic<<" "<<dynamic_left<<" "<<dynamic_l_wcet<<" "<<dynamic_h_wcet<<" "<<dynamic_l_aet<<" "<<dynamic_h_aet<<" "<<s_l;//<<std::endl;
	}
	if(expire_dynamic>0){
		
		if(dynamic_left_s>expire_dynamic){
			dynamic_left_s-=expire_dynamic;
			s_left-=expire_dynamic;
			expire_dynamic=0;
		}
		else{
			expire_dynamic-=dynamic_left_s;
			s_left-=dynamic_left_s;
			dynamic_left_s=0;
			
			if(dynamic_l_wcet_s>expire_dynamic){
				dynamic_l_wcet_s-=expire_dynamic;
				s_l-=expire_dynamic;
				expire_dynamic=0;
			}
			else{
				expire_dynamic-=dynamic_l_wcet_s;
				s_l-=dynamic_l_wcet_s;
				dynamic_l_wcet_s=0;
				if(dynamic_h_wcet_s>expire_dynamic){
					dynamic_h_wcet_s-=expire_dynamic;
					s_h-=expire_dynamic;
					expire_dynamic=0;
				}
				else{
					expire_dynamic-=dynamic_h_wcet_s;
					s_h-=dynamic_h_wcet_s;
					dynamic_h_wcet_s=0;
					if(dynamic_l_aet_s>expire_dynamic){
						dynamic_l_aet_s-=expire_dynamic;
						s_l-=expire_dynamic;
						expire_dynamic=0;
					}
					else{
						expire_dynamic-=dynamic_l_aet_s;
						s_l-=dynamic_l_aet_s;
						dynamic_h_aet_s-=expire_dynamic;
						s_h-=expire_dynamic;
						expire_dynamic=0;
					}
				}
			}
		}
		
	}
	if(expire_dynamic_total>0){
		output_slack<<"expire dynamic total "<<expire_dynamic_total<<" "<<dynamic_h_aet+static_h_aet<<"/"<<dynamic_h_wcet+static_h_wcet<<"/"<<dynamic_l_aet+static_l_aet<<"/"<<dynamic_l_wcet+static_l_wcet<<"/"<<dynamic_left+static_left<<std::endl;
		if(dynamic_left>expire_dynamic_total){
			dynamic_left-=expire_dynamic_total;
			//s_left-=expire_dynamic;
			expire_dynamic_total=0;
		}
		else{
			expire_dynamic_total-=dynamic_left;
			//s_left-=dynamic_left;
			dynamic_left=0;
			if(dynamic_l_wcet>expire_dynamic_total){
				dynamic_l_wcet-=expire_dynamic_total;
				//s_l-=expire_dynamic;
				expire_dynamic_total=0;
			}
			else{
				expire_dynamic_total-=dynamic_l_wcet;
				//s_l-=dynamic_l_wcet;
				dynamic_l_wcet=0;
				if(dynamic_h_wcet>expire_dynamic_total){
					dynamic_h_wcet-=expire_dynamic_total;
					//s_h-=expire_dynamic;
					expire_dynamic_total=0;
				}
				else{
					expire_dynamic_total-=dynamic_h_wcet;
					//s_h-=dynamic_h_wcet;
					dynamic_h_wcet=0;
					if(dynamic_l_aet>expire_dynamic_total){
						dynamic_l_aet-=expire_dynamic_total;
						//s_l-=expire_dynamic;
						expire_dynamic_total=0;
					}
					else{
						expire_dynamic_total-=dynamic_l_aet;
						//s_l-=dynamic_l_aet;
						dynamic_l_aet=0;
						dynamic_h_aet-=expire_dynamic_total;
						expire_dynamic_total=0;
					}
				}
			}

		}
		output_slack<<"expire dynamic total2 "<<expire_dynamic_total<<" "<<dynamic_h_aet+static_h_aet<<"/"<<dynamic_h_wcet+static_h_wcet<<"/"<<dynamic_l_aet+static_l_aet<<"/"<<dynamic_l_wcet+static_l_wcet<<"/"<<dynamic_left+static_left<<std::endl;
	}
	//output_slack<<std::endl;
}
void Core::handle_about_expire(){
	//deal with slack that cannot be used up by its deadline
	float dynamic_about_expire=0;
	if(slack_list.size()>0){
		std::deque<Slack>::iterator it=slack_list.begin();
		
		while(it!=slack_list.end()){
			
			if(it->expire<core_time + int(it->s/(spd[HV]-spd[LV])) ){
				//output_slack<<" | "<<it->expire<<" "<<it->s<<" "<<it->index;
				dynamic_about_expire+=it->s;
			}
			++it;
		}
	}	
	if(dynamic_about_expire>0){
		output_slack<<"re-assign slack "<<dynamic_about_expire<<"|"<<dynamic_left<<"/"<<dynamic_left_s<<"|"<<dynamic_l_wcet<<"/"<<dynamic_l_wcet_s<<"|"<<dynamic_h_wcet<<"/"<<dynamic_h_wcet_s<<"|"<<dynamic_l_aet<<"/"<<dynamic_l_aet_s<<"|"<<dynamic_h_aet<<"/"<<dynamic_h_aet_s<<std::endl;
		if(dynamic_about_expire > dynamic_left_s){
			dynamic_about_expire -= dynamic_left_s;
			dynamic_left_s += dynamic_about_expire;
			dynamic_left += dynamic_about_expire;
			s_left += dynamic_about_expire;
			if(dynamic_l_wcet_s > dynamic_about_expire){
				dynamic_l_wcet_s -= dynamic_about_expire;
				dynamic_l_wcet -= dynamic_about_expire;
				s_l -= dynamic_about_expire;
				dynamic_about_expire=0;
			}
			else{
				dynamic_about_expire -= dynamic_l_wcet_s;
				s_l -= dynamic_l_wcet_s;
				dynamic_l_wcet -= dynamic_l_wcet_s;
				dynamic_l_wcet_s=0;
				if(dynamic_h_wcet_s > dynamic_about_expire){
					dynamic_h_wcet_s -= dynamic_about_expire;
					dynamic_h_wcet -= dynamic_about_expire;
					s_h -= dynamic_about_expire;
					dynamic_about_expire=0;
				}
				else{
					dynamic_about_expire -= dynamic_h_wcet_s;
					s_h -= dynamic_h_wcet_s;
					dynamic_h_wcet -= dynamic_h_wcet_s;
					dynamic_h_wcet_s=0;
					if(dynamic_l_aet_s > dynamic_about_expire){
						dynamic_l_aet_s -= dynamic_about_expire;
						dynamic_l_aet -= dynamic_about_expire;
						s_l -= dynamic_about_expire;
						dynamic_about_expire=0;
					}
					else{
						dynamic_about_expire -= dynamic_l_aet_s;
						s_l -= dynamic_l_aet_s;
						dynamic_l_aet -= dynamic_l_aet_s;
						dynamic_l_aet_s=0;

						s_h -= dynamic_about_expire;
						dynamic_h_aet_s-=dynamic_about_expire;
						dynamic_h_aet -= dynamic_about_expire;
						dynamic_about_expire=0;
					}
				}

			}
		}
		output_slack<<"re-assign slack2 "<<dynamic_about_expire<<"|"<<dynamic_left<<"/"<<dynamic_left_s<<"|"<<dynamic_l_wcet<<"/"<<dynamic_l_wcet_s<<"|"<<dynamic_h_wcet<<"/"<<dynamic_h_wcet_s<<"|"<<dynamic_l_aet<<"/"<<dynamic_l_aet_s<<"|"<<dynamic_h_aet<<"/"<<dynamic_h_aet_s<<std::endl;
	}
}
void Core::use_slack_low(float s_used){
	float dynamic_use=s_used;
	if(dynamic_l_wcet_s>dynamic_use){
		dynamic_l_wcet_s-=dynamic_use;
		dynamic_use=0;
	}
	else{
		dynamic_use-=dynamic_l_wcet_s;
		dynamic_l_wcet_s=0;
		if(dynamic_l_aet_s>dynamic_use){
			dynamic_l_aet_s-=dynamic_use;
			dynamic_use=0;
		}
		else{
			dynamic_use-=dynamic_l_aet_s;
			dynamic_l_aet_s=0;
			static_l_total-=dynamic_use;
			dynamic_use=0;
		}
	}
	if(slack_list.size()>0){
		for(int i=0;i<slack_list.size();++i){
			if(slack_list[i].s>s_used){
				slack_list[i].s-=s_used;
				s_used=0;
				break;
			}

			else{
				s_used-=slack_list[i].s;
				slack_list[i].s=0;

			}
		}
	}
}
void Core::use_slack_high(float s_used){
	float dynamic_use=s_used;
	if(dynamic_h_wcet_s>dynamic_use){
		dynamic_h_wcet_s-=dynamic_use;
		dynamic_use=0;
	}
	else{
		dynamic_use-=dynamic_h_wcet_s;
		dynamic_h_wcet_s=0;
		if(dynamic_h_aet_s>dynamic_use){
			dynamic_h_aet_s-=dynamic_use;
			dynamic_use=0;
		}
		else{
			dynamic_use-=dynamic_h_aet_s;
			dynamic_h_aet_s=0;
			static_h_total-=dynamic_use;
			dynamic_use=0;
		}
	}
	if(slack_list.size()>0){
		for(int i=0;i<slack_list.size();++i){
			if(slack_list[i].s>s_used){
				slack_list[i].s-=s_used;
				s_used=0;
				break;
			}

			else{
				s_used-=slack_list[i].s;
				slack_list[i].s=0;

			}
		}
	}
}

void Core::use_slack_left(float s_used){
	float dynamic_use=s_used;
	if(dynamic_left_s>dynamic_use){
		dynamic_left_s-=dynamic_use;
		dynamic_use=0;
	}
	else{
		dynamic_use-=dynamic_left_s;
		dynamic_left_s=0;
		static_left_total-=dynamic_use;
		dynamic_use=0;
	}
	if(slack_list.size()>0){
		for(int i=0;i<slack_list.size();++i){
			if(slack_list[i].s>s_used){
				slack_list[i].s-=s_used;
				s_used=0;
				break;
			}

			else{
				s_used-=slack_list[i].s;
				slack_list[i].s=0;

			}
		}
	}
}

void Core::set_dvfs_mode(bool dvfs_type){
	dvfs_tradition=dvfs_type;

}

Task* copy_task(Task base, int j , int actual){
	Task* temp_task=new Task;
	temp_task->name=base.name;
	temp_task->index=base.index;
	temp_task->period=base.period;
	temp_task->wcet=base.wcet;
	temp_task->arrival=base.arrival+j*base.period;
	temp_task->actual=actual;
	temp_task->deadline=base.deadline+j*base.period;
	temp_task->exe=base.exe;
	temp_task->mean=base.mean;
	temp_task->high=base.high;
	temp_task->read_ratio=base.read_ratio;
	return temp_task;
}
long get_gcd(long a, long b){
	if(b==0){
		return a;
	}
	else{
		return get_gcd(b, a%b);
	}
}
long get_lcm(long a, long b){
	//std::cout<<"gcd is"<<get_gcd(a,b)<<std::endl;
	return a*b/get_gcd(a,b);
}
bool cmp_deadline(const Task& a, const Task& b){
	return a.deadline<b.deadline;
}
bool cmp_arrival(const Task& a, const Task& b){
	return a.arrival<b.arrival;
}
bool cmp_slack(const Slack& a, const Slack& b){
	return a.expire<b.expire;
}
bool defer_cmp_deadline_reverse(const Defer_task& a, const Defer_task& b){
	return a.deadline>b.deadline;
}
/*Task* temp_task=new Task;
			temp_task->name=task_list[0].name;
			temp_task->index=task_list[0].index;
			temp_task->period=task_list[0].period;
			temp_task->wcet=task_list[0].wcet;
			temp_task->arrival=task_list[0].arrival;
			temp_task->actual=task_list[0].actual;
			temp_task->deadline=task_list[0].deadline;
			temp_task->exe=0;
			temp_task->mean=task_list[0].mean;*/
/*Task* temp_task=new Task;
			temp_task->name=t_list[i].name;
			temp_task->index=t_list[i].index;
			temp_task->period=t_list[i].period;
			temp_task->wcet=t_list[i].wcet;
			temp_task->arrival=t_list[i].arrival+j*t_list[i].period;
			temp_task->actual=t_list[i].wcet/10*8;
			temp_task->deadline=t_list[i].deadline+j*t_list[i].period;
			temp_task->exe=task_list[i].exe;*/
/*int count=0;
	while (count<5){
		getline(*files[task_list[0].index],power_string);
		std::cout<<task_list[0].index<<" "<<power_string<<std::endl;
		count+=1;
	}
	count=0;
	while (count<5){
		getline(*files[task_list[1].index],power_string);
		std::cout<<task_list[1].index<<" "<<power_string<<std::endl;
		count+=1;
	}
	count=0;
	while (count<5){
		getline(*files[task_list[0].index],power_string);
		std::cout<<task_list[0].index<<" "<<power_string<<std::endl;
		count+=1;
	}
	files[task_list[0].index]->close();
	files[task_list[0].index]->open(task_list[0].name.c_str());
	count=0;
	while (count<5){
		getline(*files[task_list[0].index],power_string);
		std::cout<<task_list[0].index<<" "<<power_string<<std::endl;
		count+=1;
	}*/
	
	/*
	int count=0;
	for (int i=0;i<100;i++){
		if(count<10){
			getline(file,power_string);
			std::istringstream line(power_string);
			line>>IPC>>power[0]>>power[1]>>power[2]>>power[3]>>power[4]>>power[5]>>power[6]>>power[7]>>power[8]>>power[9]>>power[10]>>power[11]>>power[12]>>power[13]>>power[14];
			for(int i=0;i<blocks;++i){
				power_cal[i]+=power[i];
			}
			count+=1;
		}
		else{
			for(int i=0;i<blocks;++i){
				power_cal[i]/=10.0;
			}
			compute_temp_tilts(tmpr,power_cal,TIME_SLICE,N_BLOCKS,N_BLOCKS*NL+EXTRA,Aj,Bx);
			std::cout<<power_cal[3]<<" "<<tmpr[3]<<std::endl;
			for(int i=0;i<EXTRA+blocks;++i){
				power_cal[i]=0;
			}
			count=0;
			
		}
	}*/