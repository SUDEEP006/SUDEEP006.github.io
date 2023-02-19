clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%############# SIGNAL AQUISITION (FORWARD PROBLEM) ###################

%%%%%%%% 2D Computational grid %%%%%%%%%%%%
Nx = 129; Ny = 129;
dx = 1e-3; dy = 1e-3;
kx = 65; ky = 65;
ROIx=40; ROIy=40;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%%%%%%%% Medium %%%%%%%%%%%%
speed_of_sound = 1500; % [m/s]
density11 = 1000; % [Kg/m^3]
medium.density = density11;
medium.sound_speed = speed_of_sound*ones(Nx, Ny);

%%%%%%%%%%%%%% Time array  %%%%%%%%%%%%%%
[kgrid.t_array dt] = makeTime(kgrid,medium.sound_speed);
time_1 = kgrid.t_array;

%%%%%%%% Source %%%%%%%%%%%%
% % disc1_x_pos = kx;  % [grid points]
% % disc1_y_pos = ky;  % [grid points]
% % disc1_radius = 10; % [grid points]
% % press_mag1 = 1;     % [Pa]
% % disc1 = press_mag1*makeDisc(Nx, Ny, disc1_x_pos, disc1_y_pos, disc1_radius);

img = loadImage('vasculature.png');
img = imresize(img,[Nx Ny]);
source.p0 = img;

%%%%%%%% Sensor %%%%%%%%%%%%
sensor_radius = 60e-3;
num_sensor_points = 100;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

center_freq = 2.25e6;                                  % [Hz] Frequency response
bandwidth = 70;

signal_to_noise_ratio = 40;                            % [dB] Noise

PML_size = [20,20];                                    % Perfect Matching Layer (PML)
PML_Alpha = [10,10];
input_args = {'PMLSize',PML_size,'PMLAlpha',PML_Alpha,'PMLInside',false,'PlotPML', true}; %input arguments

%%%%%%%%%%%%%%define the frequency response of the sensor elements
sensor.frequency_response = [center_freq, bandwidth];

%%%%%%%%%%%%%% Run the simulation %%%%%%%%%%%%%%
sensor_data_filtered = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%%%%%%%%%%%%%%add noise to the recorded sensor data
sensor_data = addNoise(sensor_data_filtered, signal_to_noise_ratio, 'peak');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%############# RECONSTRUCTION (INVERSE PROBLEM) ###################

%%%%%%%%%%%%%% back projection formula %%%%%%%%%%%%%%
nn = length(sensor_data(1,:));
del_p = zeros(1,nn);
b = zeros(num_sensor_points, nn);
for n = 1:1:num_sensor_points
    p_sensor = sensor_data(n,:);
    %     diff_p_sensor = diff(p_sensor);
    %     del_p(1:nn-1) = diff_p_sensor;
    %     b(n,:) = 2 * p_sensor - (2 * time_1 .* del_p)/dt;
    b(n,:) = 2 * p_sensor;
end

%%%%%%%%%%%%%%% Center Coordinates %%%%%%%%%%%%%%
x_center = kx*dx;
y_center = ky*dy;

delta_theta = 360/num_sensor_points;                       % Angular positions
theta = (-180:delta_theta:179)*pi/180;

%%%%%%%%%%%%%% Sensor Coordinates %%%%%%%%%%%%%%
[x_sensor] = sensor_radius*sin(theta) + x_center;
[y_sensor] = sensor_radius*cos(theta) + y_center;


%%%%%%%%%%%%%% A grid of zeroes having same size as computational grid
image_recon = zeros(Nx, Ny);

for i11= kx-ROIx:1:kx+ROIx
    for j11=ky-ROIy:1:ky+ROIy
        for n=1:1:num_sensor_points
            distance11 = sqrt((x_sensor(n)- i11*dx).^2 + (y_sensor(n)- j11*dy).^2);
            ntcell = floor(1 + distance11/(speed_of_sound * dt));
            image_recon(i11,j11) = image_recon(i11,j11) + b(n,ntcell)/(distance11*distance11);
        end
    end
end

reconstructed_image = image_recon;
reconstructed_image = reconstructed_image/max(max(reconstructed_image));

figure(1)
imagesc(source.p0);
colormap(gray)
colormap(flipud(gray))
colorbar('westoutside')
axis image;

figure(2)
imagesc(reconstructed_image);
colormap(gray)
colormap(flipud(gray))
colorbar('westoutside')
axis image;

