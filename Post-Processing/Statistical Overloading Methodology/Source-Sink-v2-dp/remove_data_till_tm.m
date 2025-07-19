function [id,dp,tm]...
    = remove_data_till_tm(id,dp,tm,time_gap)

cond = tm>=time_gap;
id = id(cond);
dp = dp(cond);
tm = tm(cond);