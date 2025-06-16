%> @file load_cst_demo.m
%> @brief Visualizes CST far-field simulation data in MATLAB using custom coordinate conversion.
%>
%> This script demonstrates how to load and visualize far-field results exported from CST Studio Suite.
%> It compares the raw CST spherical data to the transformed version in the simulation’s custom coordinate
%> system using the `DishAnalyzer` plotting tools.
%>
%> @section inputs Inputs
%> - CST far-field export file: 'farfield_(f=2.4)_[1].txt'
%> - Function `loadcstfarfield`: loads raw CST spherical coordinate data
%> - Function `loadcstfarfield_mycoords`: loads CST data and converts it to custom spherical coordinates
%>
%> @section outputs Outputs
%> - 3D surface plots of |E| in dB for both raw CST and custom coordinate systems
%>
%> @section usage Usage
%> - Requires `DishAnalyzer` class for 3D plotting
%> - Designed for quick inspection of CST far-field formats and transformations
%>
%> @note This script does not perform analysis—only visualization. Data fidelity must be ensured upstream.


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