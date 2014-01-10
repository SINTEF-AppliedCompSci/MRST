do_plot = true;
hold on
figure;
[T_w, T_o] = readUpscaledRelperm(do_plot, 'fixed', false); hold on

[T_w, T_o] = readUpscaledRelperm(do_plot, 'fixed', true);
%figure
[T_w, T_o] = readUpscaledRelperm(do_plot, 'bc', false);
[T_w, T_o] = readUpscaledRelperm(do_plot, 'bc', true);



