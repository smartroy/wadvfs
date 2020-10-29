#ifndef CORE_H
#define CORE_H
#include <vector>
#include <deque>
#include <string>
#include <fstream>
#include <sstream>

extern "C"{
	#include "flp.h"
	#include "temperature.h"
	#include "tilts.h"
}
const int blocks(15);
const int dvfs_step(10);
const int LV(0);
const int MV(1);
const int HV(2);
const double TIME_STEP(5e-3);
//const int EXTRA(12);
//const int NL(4);
struct Task{
	std::string name;
	short index;
	int period;
	int wcet;
	long arrival;
	int actual;
	long deadline;
	int exe; 
	float mean;
	float high;
	float read_ratio;
	//std::ifstream file;
};
struct Slack{
	float s;
	int expire;
	int index;
};
struct Defer_task{
	float left;
	int deadline;
	int index;
	float utilization;
};
class Core{
private:
	
	std::deque<Task> exe_list;
	std::deque<Task> task_list;
	std::deque<Slack> slack_list;
	std::deque<Defer_task> defer_list;
	long core_lcm;
	long sys_lcm;
	flp_t *flp;
	RC_model_t *rc_model;
	thermal_config_t config;
	double Aj[24576];
	double Bx[24576];
	float s_l,s_h,s_left;
	float s_h_aet,s_l_aet;
	float s_h_wcet,s_l_wcet;
	float static_l_aet,static_h_aet,dynamic_l_aet,dynamic_h_aet;
	float static_l_wcet,static_h_wcet,dynamic_l_wcet,dynamic_h_wcet;
	float dynamic_h_aet_s,dynamic_l_aet_s,dynamic_h_wcet_s,dynamic_l_wcet_s,dynamic_left_s;
	float static_h_total,static_l_total,static_left_total;
	float s_static,s_dynamic,static_total;
	float static_left;
	float dynamic_left;
	float core_utilization;	
	float wc_utilization;

	float exe_h,exe_l;
	float s_static_l,s_static_h;
	float s_dynamic_h,s_dynamic_l;

	int spd_h;
	int spd_l;
	int voltage;
	float core_high;
	int core_wcet;
	float core_low;
	float new_dynamic;
	int spd[HV+1];
	double p_ratio[HV+1];
	double area_ratio[blocks];
	bool dvfs_tradition;
	float task_utilization[10];
	struct Defer_task defer_array[10];
	float dynamic_array[10];
	float static_array[10];
	float total_power;
	int task_count;
	

public:
	double temperature[blocks];
	double power[blocks];
	double power_cal[blocks+EXTRA];
	double tmpr[NL*blocks+EXTRA];
	double IPC;
	double IPC_dvfs;
	long core_time;
	int core_speed;
	int test_power;
	float ipc_threshold;
	double rel_em;
	double rel_bd;
	double rel_cycle;
	double rel_total;
	float expired_slack_total;

	Core();
	void set_temperature();
	void set_temperature(float *temp_array);
	void dis_temperature();
	void gen_exe();
	void add_task(Task new_task);
	void gen_lcm(std::vector<Task> t_list,float threshold);
	void gen_flp();
	void print_RC();
	void run_lcm(long sys_wide_lcm,double sys_time);
	void core_read(int read_speed);
	void dynamic_reset();
	void cal_static();
	void alloc_dynamic();
	void handle_expire();
	void cal_rel(double sys_time);
	void set_dvfs_mode(bool dvfs_type);
	void handle_expire_tradition();
	void alloc_dynamic_tradition();
	void cal_static_tradition();
	void set_speed();
	void pdvfs_defer();
	void defer_set_speed(float ratio);
	void use_slack_high(float s_used);
	void use_slack_low(float s_used);
	void use_slack_left(float s_used);
	void handle_about_expire();
};
long get_gcd(long a, long b);
long get_lcm(long a, long b);
bool cmp_deadline(const Task& a, const Task& b);
bool cmp_arrival(const Task& a, const Task& b);
Task* copy_task(Task base, int j,int actual );
bool cmp_slack(const Slack& a, const Slack& b);
bool defer_cmp_deadline_reverse(const Defer_task& a, const Defer_task& b);
#endif