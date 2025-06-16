[EdB_, Etheta_, Ephi_, THETA_, PHI_] = loadcstfarfield('farfield_(f=2.4)_[1].txt');

dish_analyzer = DishAnalyzer([]);

figure;
dish_analyzer.plot_3d_rad_pattern([], EdB_, THETA_, PHI_);
title("Raw data from CST");
% set(gcf,'Position',[100 100 600 600])

%%
[EdB_, Etheta_, Ephi_, THETA_, PHI_] = loadcstfarfield_mycoords('farfield_(f=2.4)_[1].txt');

dish_analyzer = DishAnalyzer([]);

figure;
dish_analyzer.plot_3d_rad_pattern([], EdB_, THETA_, PHI_);
title("Data from CST, Converted to Our Coordinate System");



% figure;
% surf(rad2deg(THETA), rad2deg(PHI), EdB, 'EdgeColor', 'none');
% xlabel('\theta (deg)'); ylabel('\phi (deg)'); title('|E_\theta| (dB)');
% colorbar;