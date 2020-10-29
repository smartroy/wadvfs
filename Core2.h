#ifndef CORE2_H
#define CORE2_H
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <sstream>
extern "C"
{
#include "flp.h"
#include "temperature.h"
#include "tilts.h"
}

const double SCALE_EM(0.3);
const double SCALE_BD(0.18);
const double K_BOL(8.62e-5);
const double EA_EM(0.9);
const double EXP_EM(2.0);
const double EXP_BD(2.0);

const int blocks(15);
const int dvfs_step(10);
const int LV(0);
const int MV(1);
const int HV(2);
const double TIME_STEP(5e-3);
struct BaseTask
{
    std::string name;
    int index;
    int period;
    int wcet;
    float mean;
    float high;
    
    float utilizaton;
    //std::ifstream file;
};

struct Task
{
    long arrival;
    int actual;
    long deadline;
    int exe;
    float utilization;
    int index;
    float read_ratio;
    int wcet;
    Task(int offset, const BaseTask &t);
};

enum SlackType
{
    eStatic,
    eDynamic
};

struct Slack
{
    float s;
    int expire;
    int index;
    SlackType type;
};

struct CompareTaskDeadline
{
    bool operator()(const Task &t1, const Task &t2)
    {
        return t1.deadline > t2.deadline;
    }
};

struct CompareTaskArrival
{
    bool operator()(const Task &t1, const Task &t2)
    {
        return t1.arrival > t2.arrival;
    }
};

struct CompareSlackExipre
{
    bool operator()(const Slack &s1, const Slack &s2)
    {
        if (s1.expire == s2.expire)
        {
            return s1.type < s2.type;
        }
        return s1.expire > s2.expire;
    }
};

class Core
{

    std::priority_queue<Task, std::vector<Task>, CompareTaskDeadline> exe_list;
    std::priority_queue<Task, std::vector<Task>, CompareTaskArrival> task_list;
    std::priority_queue<Slack, std::vector<Slack>, CompareSlackExipre> slack_list;
    std::vector<BaseTask> assigned_tasks;
    long core_lcm;
    long sys_lcm;
    flp_t *flp;
    RC_model_t *rc_model;
    thermal_config_t config;
    double Aj[24576];
    double Bx[24576];

    float slackLow, slackHigh, slackFree; // total slack available for low/high IPC and slack can be used for both
    float core_utilization;               //current core utlization, for pdvfs
    float wc_utilization;                 // worst case utilization

    float highIpcPhase; // high/low IPC execution time if run to WCET
    float lowIpcPhase;
    float highIpcPhaseAet;
    float lowIpcPhaseAet;
    float highIpcPhaseWcet;
    float lowIpcPhaseWcet;

    float slackAssignHigh, slackAssignLow; //assigned slack for low/high
    float slackNeedHigh, slackNeedlow;
    int voltage;
    int spd[HV + 1];
    double p_ratio[HV + 1];    //power ratio when using different voltage
    double area_ratio[blocks]; // weight for each functional block when calculating reliability
    bool dvfs_tradition;
    bool new_release;
    float total_power;
    int task_count;
    std::vector<BaseTask> coreTask;
    float task_utilization[10];

public:
    double temperature[blocks];
    double power[blocks];
    double power_cal[blocks + EXTRA];
    double tmpr[NL * blocks + EXTRA];
    double IPC;
    double IPC_dvfs; //IPC value after a decision step, used for determine if dvfs
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
    void prep_core(std::vector<BaseTask>, float);
    void gen_lcm();
    void reset_slack();
    void gen_flp();
    void print_RC();
    void run_lcm(long sys_wide_lcm, double sys_time);
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
#endif