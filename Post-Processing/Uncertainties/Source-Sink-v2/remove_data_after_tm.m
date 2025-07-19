function [id,dp,tm] = remove_data_after_tm(id,dp,tm,time)

idx = (tm>time);
tm(idx) = [];
dp(idx) = [];
id(idx) = [];