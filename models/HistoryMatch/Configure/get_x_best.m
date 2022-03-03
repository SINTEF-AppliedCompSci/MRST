function [x_best, fval_sum_best] = get_x_best(model)

x_list = model.history_match.x_list;
fval_sum = sum(model.history_match.fval_list,2);
sorted_fval_sum = sort(fval_sum);
fval_sum_best = fval_sum(fval_sum==sorted_fval_sum(1),:);
x_best = x_list(fval_sum==sorted_fval_sum(1),:);