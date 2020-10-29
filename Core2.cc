#include "Core2.h"
#include <iostream>

#include <cmath>
#include <iomanip>

#include <algorithm>
#include <chrono>
#include <random>
#include <assert.h>

extern std::vector<std::ifstream *> files;
extern std::ofstream output;
extern std::ofstream output_slack;
extern std::ofstream output_tmpr;
extern std::ofstream output_rel;
extern std::ofstream output_ipc;
extern std::ofstream output_dvfs;

Core::Core()
{
    set_temperature();
    gen_flp();
}

void Core::set_temperature()
{
    for (int i = 0; i < blocks; i++)
    {
        temperature[i] = T_INIT + 273.15;
        power[i] = 0.0;
        power_cal[i] = 0.0;
        tmpr[i] = T_INIT + 273.15;
        area_ratio[i] = 1;
    }
    for (int i = blocks; i < EXTRA + blocks; i++)
    {
        power_cal[i] = 0.0;
    }
    for (int i = blocks; i < NL * blocks + EXTRA; i++)
    {
        tmpr[i] = T_INIT + 273.15;
    }

    spd[LV] = 6;
    spd[MV] = 6;
    spd[HV] = 10;
    p_ratio[LV] = 0.41; //0.41;//0.32;0.54;0.68;
    p_ratio[MV] = 0.41;
    p_ratio[HV] = 1.0;
    rel_em = 1;
    rel_bd = 1;
    rel_cycle = 1;
    rel_total = 1;
    expired_slack_total = 0;
}

void Core::set_temperature(float *temp_array)
{
    for (int i = 0; i < blocks; ++i)
    {
        temperature[i] = temp_array[i];
    }
}

void Core::dis_temperature()
{
    for (int i = 0; i < blocks; i++)
    {
        std::cout << temperature[i] << std::endl;
    }
}

void Core::gen_flp()
{
    flp = read_flp("90nm.flp", FALSE);
    config = default_thermal_config();
    rc_model = alloc_RC_model(&config, flp);
    populate_R_model(rc_model, flp);
    populate_C_model(rc_model, flp);
    initialize_tilts(rc_model, flp, rc_model->config->sampling_intvl, Aj, Bx);
}
void Core::print_RC()
{
    debug_print_model(rc_model);
}

void Core::reset_slack()
{
    slackLow = 0;
    slackHigh = 0;
    slackFree = 0;

    highIpcPhase = 0; // high/low IPC execution time if run to WCET
    lowIpcPhase = 0;
    highIpcPhaseAet = 0;
    lowIpcPhaseAet = 0;
    highIpcPhaseWcet = 0;
    lowIpcPhaseWcet = 0;

    slackAssignHigh = 0;
    slackAssignLow = 0;
    slackNeedHigh = 0;
    slackNeedlow = 0;
}

void Core::prep_core(std::vector<BaseTask> tasks, float threshold)
{
    assigned_tasks = tasks;
    ipc_threshold = threshold;
}

void Core::gen_lcm()
{
    assert(assigned_tasks.size() > 0);

    wc_utilization = 0;

    reset_slack();

    core_lcm = assigned_tasks[0].period;

    for (int i = 0; i < assigned_tasks.size(); i++)
    {
        core_lcm = get_lcm(core_lcm, assigned_tasks[i].period);
        wc_utilization += assigned_tasks[i].utilizaton;
    }
    core_utilization = wc_utilization;

    while (!task_list.empty())
        task_list.pop();
    for (int i = 0; i < assigned_tasks.size(); i++)
    {
        for (int j = 0; j < core_lcm / assigned_tasks[i].period; j++)
        {
            task_list.push(Task(j, assigned_tasks[i]));
        }
    }
}

void Core::run_lcm(long sys_wide_lcm, double sys_time)
{
    sys_lcm = sys_wide_lcm;
    core_time = 0;
    core_speed = spd[HV];
    unsigned seed = 5;
    std::default_random_engine generator(seed);
    bool new_release = false;
    float out_ipc = 0;

    while (core_time < sys_lcm)
    {
        if (core_time % core_lcm == 0)
        {
            gen_lcm();
        }
        while (!task_list.empty() && task_list.top().arrival <= core_time)
        {
            std::normal_distribution<double> distibution(assigned_tasks[task_list.top().index].mean, 0.2);
        }
    }
}

long get_gcd(long a, long b)
{
    if (b == 0)
    {
        return a;
    }
    else
    {
        return get_gcd(b, a % b);
    }
}
long get_lcm(long a, long b)
{
    //std::cout<<"gcd is"<<get_gcd(a,b)<<std::endl;
    return a * b / get_gcd(a, b);
}