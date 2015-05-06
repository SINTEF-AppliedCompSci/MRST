function addLegends(ax, leg_linecolors, linecolors, leg_linetypes, linetypes)
   h1 = ax;
   h3 = axes('Position', get(h1, 'Position'));
   x = 1:2;

   % linw colors
   assert(numel(leg_linecolors) == numel(linecolors))
   leg_names = {leg_linecolors{1:numel(linecolors)}};%#ok
   tmp = nan(numel(x), numel(leg_names));
   f2 = plot(x, tmp); % dummy plot
   axis off;
   set(gca, 'XtickLabel', [])

   for f = 1:numel(f2)
      set(f2(f), 'Marker', 'none');
      set(f2(f), 'MarkerFaceColor', 'none');
      set(f2(f), 'LineStyle', ' - '); % markers{f});
      set(f2(f), 'Color', linecolors{f});
   end

   set(h3, 'FontSize', 16);
   bb = legend(leg_names, 'Location', 'SouthEast', 'FontSize', 16);
   loc = get(bb, 'Position');
   loc(2) = loc(2) + 0.15;
   set(bb, 'Position', loc)
   set(gca, 'FontSize', 19)


   h4 = axes('Position', get(h1, 'Position'));%#ok

   assert(numel(leg_linetypes) == numel(linetypes));
   leg_names = {leg_linetypes{1:numel(linetypes)}};%#ok

   x = 1:2;
   tmp = nan(numel(x), numel(leg_names));
   f2 = plot(x, tmp);
   axis off;
   set(gca, 'XtickLabel', []);

   for f = 1:numel(f2)
      set(f2(f), 'Marker', 'none');
      set(f2(f), 'MarkerFaceColor', 'none');
      set(f2(f), 'LineStyle', linetypes{f});
      set(f2(f), 'Color', 'k');
   end

   set(h3, 'FontSize', 16);
   bb = legend(leg_names, 'Location', 'SouthEast', 'FontSize', 16);
   loc = get(bb, 'Position');%#ok
   set(gca, 'FontSize', 19)

   set(gcf,'PaperPositionMode','auto');
end
