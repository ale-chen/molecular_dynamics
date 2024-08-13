% Simulation parameters
num_bodies = 3;
total_time = 2;
dt = 0.001;
k_inter = 80;
k_center = 3;
sigma_micro = 2;

% Create rigid bodies with randomized initial positions
rigid_bodies = cell(num_bodies, 1);
for i = 1:num_bodies
    disp(['Drawing Rigid Body ', num2str(i), ' of ', num2str(num_bodies), ':']);
    rigid_bodies{i} = create_rigid_body();
    close;  % Close the figure after each rigid body is drawn
    disp(['Rigid Body ', num2str(i), ' created successfully.']);
end

% Ensure non-overlapping initial positions
rigid_bodies = initialize_positions(rigid_bodies, sigma_micro);

% Initialize velocities
for i = 1:num_bodies
    rigid_bodies{i}.velocity_history(1, :) = (rand(1, 2) - 0.5) * 0.1;
    rigid_bodies{i}.angular_velocity_history(1) = (rand(1) - 0.5) * 0.1;
end

% Initialize time
t = 0;

% Calculate total number of simulation steps
total_steps = round(total_time / dt);

% Initialize progress bar
disp('Starting simulation:');
progress_bar_width = 50;
fprintf(['[' repmat(' ', 1, progress_bar_width) ']\n\n']);

% Simulation loop
for step = 1:total_steps
    % Update positions and orientations of rigid bodies
    for i = 1:num_bodies
        rigid_bodies{i} = update_rigid_body(rigid_bodies{i}, dt);
    end
    
    % Calculate forces and torques on each rigid body
    for i = 1:num_bodies
        [force_total, torque_total] = calculate_total_force_torque(rigid_bodies, i, k_inter, k_center, sigma_micro);
        rigid_bodies{i} = apply_force_torque(rigid_bodies{i}, force_total, torque_total);
    end
    
    % Increment time
    t = t + dt;
    
    % Update progress bar
    progress = step / total_steps;
    filled_width = round(progress * progress_bar_width);
    fprintf('\033[1A');  % Move cursor up one line
    fprintf('\033[K');   % Clear the line
    fprintf(['[' repmat('=', 1, filled_width) repmat(' ', 1, progress_bar_width - filled_width) '] %.1f%%\n'], progress * 100);
    
    % Print checkpoint every 10% of simulation
    if mod(step, round(total_steps / 10)) == 0
        disp(['Simulation ', num2str(progress * 100), '% complete. Current time: ', num2str(t), ' seconds.']);
    end
end

disp('Simulation completed. Starting visualization...');

% Visualization
figure;
hold on;
axis equal;
xlim([-50, 50]);
ylim([-50, 50]);
title('Rigid Body Simulation');

num_frames = size(rigid_bodies{1}.position_history, 1);

% Initialize progress bar for visualization
disp('Visualizing simulation:');
fprintf(['[' repmat(' ', 1, progress_bar_width) ']\n\n']);

for frame = 1:num_frames
    clf;
    hold on;
    for i = 1:num_bodies
        coords = get_particle_coordinates(rigid_bodies{i}, frame);
        plot(coords(:, 1), coords(:, 2), 'o');
    end
    xlim([-50, 50]);
    ylim([-50, 50]);
    title(['Time: ', num2str((frame-1)*dt)]);
    drawnow;
    pause(0.01);
    
    % Update progress bar
    progress = frame / num_frames;
    filled_width = round(progress * progress_bar_width);
    fprintf('\033[1A');  % Move cursor up one line
    fprintf('\033[K');   % Clear the line
    fprintf(['[' repmat('=', 1, filled_width) repmat(' ', 1, progress_bar_width - filled_width) '] %.1f%%\n'], progress * 100);
    
    % Print checkpoint every 10% of visualization
    if mod(frame, round(num_frames / 10)) == 0
        disp(['Visualization ', num2str(progress * 100), '% complete. Current frame: ', num2str(frame), ' of ', num2str(num_frames)]);
    end
end

disp('Visualization completed.');
% Helper functions

function rigid_body = create_rigid_body()
    disp('Please draw a polygon.');
    disp('Left-click to add points, right-click or press Enter to finish drawing.');
    figure;
    hold on;
    axis equal;
    xlim([-10 10]);
    ylim([-10 10]);
    sigma = 2;

    x = [];
    y = [];
    while true
        [xi, yi, button] = ginput(1);
        if button == 1  % Left-click
            x = [x; xi];
            y = [y; yi];
            plot(x, y, 'b-', 'LineWidth', 2);
        else  % Right-click, Enter key, or any other button
            break;
        end
    end

    [x_opt, y_opt] = optimize_vertices(x, y, sigma);
    particle_centers = generate_particle_centers(x_opt, y_opt, sigma);

    mass_vector = 0.5 + rand(size(particle_centers, 1), 1);
    rigid_body = RigidBody(particle_centers, mass_vector);
end

function [x_opt, y_opt] = optimize_vertices(x, y, sigma)
    num_vertices = length(x);
    
    % Define objective
    fun = @(coords) vertex_objective(coords, x, y);
    
    % Define constraint
    nonlcon = @(coords) vertex_constraints(coords, sigma, num_vertices);
    
    % Set initial guess
    coords_init = [x; y];
    
    % Set optimization options
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
    
    % Perform constrained optimization
    coords_opt = fmincon(fun, coords_init, [], [], [], [], [], [], nonlcon, options);
    
    % Extract optimized x and y coordinates
    x_opt = coords_opt(1:num_vertices);
    y_opt = coords_opt(num_vertices+1:end);
end

function obj = vertex_objective(coords, x, y)
    num_vertices = length(x);
    x_adj = coords(1:num_vertices);
    y_adj = coords(num_vertices+1:end);
    
    % Calculate sum of squared distances between original and adjusted vertices
    obj = sum((x_adj - x).^2 + (y_adj - y).^2);
end

function [c, ceq] = vertex_constraints(coords, sigma, num_vertices)
    x_adj = coords(1:num_vertices);
    y_adj = coords(num_vertices+1:end);
    
    % Calculate edge lengths
    edge_lengths = sqrt((x_adj(2:end) - x_adj(1:end-1)).^2 + (y_adj(2:end) - y_adj(1:end-1)).^2);
    
    % Inequality constraints (edge lengths >= sigma)
    c = sigma - edge_lengths;
    
    % Equality constraints (edge lengths == integer multiples of sigma)
    ceq = mod(edge_lengths, sigma);
end

function centers = generate_particle_centers(x, y, sigma)
    num_vertices = length(x);
    centers = [x, y];
    
    for i = 1:num_vertices
        j = mod(i, num_vertices) + 1;
        edge_length = sqrt((x(j) - x(i))^2 + (y(j) - y(i))^2);
        num_particles_edge = round(edge_length / sigma) - 1;
        
        if num_particles_edge > 0
            t = linspace(0, 1, num_particles_edge+2);
            edge_coords = [linspace(x(i), x(j), num_particles_edge+2)', ...
                           linspace(y(i), y(j), num_particles_edge+2)'];
            centers = [centers; edge_coords(2:end-1, :)];
        end
    end
end

function rigid_bodies = initialize_positions(rigid_bodies, sigma_micro)
    num_bodies = length(rigid_bodies);
    min_distance = 10 * sigma_micro;  % Adjust the minimum distance based on macro particle size
    
    for i = 1:num_bodies
        while true
            % Generate random initial position
            pos = (rand(2, 1) - 0.5) * 100;
            rigid_bodies{i}.position_history(1, :) = pos';
            
            % Check for overlap with other macro particles
            overlap = false;
            for j = 1:i-1
                dist = norm(pos - rigid_bodies{j}.position_history(1, :)');
                if dist < min_distance
                    overlap = true;
                    break;
                end
            end
            
            if ~overlap
                break;
            end
        end
    end
end

function rigid_body = update_rigid_body(rigid_body, dt)
    % Update position and orientation
    rigid_body.position_history(end+1, :) = rigid_body.position_history(end, :) + rigid_body.velocity_history(end, :) * dt;
    rigid_body.orientation_history(end+1) = rigid_body.orientation_history(end) + rigid_body.angular_velocity_history(end) * dt;
end

function [force_total, torque_total] = calculate_total_force_torque(rigid_bodies, body_index, k_inter, k_center, sigma_micro)
    num_bodies = length(rigid_bodies);
    force_total = zeros(2, 1);
    torque_total = 0;
    
    % Interaction forces between particles of different rigid bodies
    for j = 1:num_bodies
        if body_index ~= j
            coords_i = get_particle_coordinates(rigid_bodies{body_index}, size(rigid_bodies{body_index}.position_history, 1));
            coords_j = get_particle_coordinates(rigid_bodies{j}, size(rigid_bodies{j}.position_history, 1));
            
            for k = 1:size(coords_i, 1)
                for l = 1:size(coords_j, 1)
                    diff = coords_i(k, :)' - coords_j(l, :)';
                    distance = norm(diff);
                    
                    if distance <= sigma_micro
                        direction = diff / distance;
                        force = k_inter * direction;
                        
                        force_total = force_total + force;
                        r = coords_i(k, :) - rigid_bodies{body_index}.center_of_mass;
                        torque = cross([r, 0], [force', 0]);  % Ensure 3D vectors for cross product
                        torque_total = torque_total + torque(3);
                    end
                end
            end
        end
    end
    
    % Central potential force
    force_center = -k_center * rigid_bodies{body_index}.position_history(end, :)';
    force_total = force_total + force_center;
end

function rigid_body = apply_force_torque(rigid_body, force_total, torque_total)
    % Apply total force and torque to the rigid body
    rigid_body.velocity_history(end+1, :) = rigid_body.velocity_history(end, :) + (force_total' / sum(rigid_body.mass_vector));
    rigid_body.angular_velocity_history(end+1) = rigid_body.angular_velocity_history(end) + torque_total / rigid_body.moment_of_inertia;
end

function coords = get_particle_coordinates(rigid_body, frame_index)
    rotation_matrix = [cos(rigid_body.orientation_history(frame_index)), -sin(rigid_body.orientation_history(frame_index)); 
                       sin(rigid_body.orientation_history(frame_index)), cos(rigid_body.orientation_history(frame_index))];
    relative_positions = rigid_body.particle_centers - repmat(rigid_body.center_of_mass, size(rigid_body.particle_centers, 1), 1);
    rotated_positions = (rotation_matrix * relative_positions')';
    coords = rotated_positions + repmat(rigid_body.position_history(frame_index, :), size(rigid_body.particle_centers, 1), 1);
end

function moment_of_inertia = calculate_moment_of_inertia(centers, mass_vector)
    moment_of_inertia = sum(mass_vector .* sum(centers.^2, 2));
end