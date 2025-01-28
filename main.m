% we rewrite tessellation based on topology and rotation bwtween initial
% pattern and deployed pattern.

clearvars
addpath(genpath('.'))
shape_name = 'circle';

% Set the type of the initial map for the optimization in the deployed space
% 1: standard deployed configuration
% 2: standard deployed configuration with rescaling (with optional parameter scale_factor)
% 3: conformal map (Schwarz-Christoffel mapping)
% 4: Teichmuller map (Meng et al., SIIMS 2016)
% initial_map_type = 1;
% initial_map_type = 2; scale_factor = 0.75; % or other positive number
initial_map_type = 3;
% initial_map_type = 4;

if ~exist('scale_factor','var')
    scale_factor = [];
end
% Require the contracted pattern to be rectangular? 0 (no) / 1 (yes)
% fix_contracted_boundary_shape = 0;
fix_contracted_boundary_shape = 1;

% Further specify the width-to-height ratio if the pattern is required to 
% be rectangular? 0 (no/not applicable) / other positive number (the prescribed ratio)
% rectangular_ratio = 0; % no need to specify / not applicable
rectangular_ratio = 1; % desired to be a square
% rectangular_ratio = 2; % 2 means desired width-to-height ratio = 2:1, can be changed to other values

width = 4;
height = 4;
plot_toggle = true;
[pointsD_standard, face_setsD, edgesD, edge_pairsD, anglesD, int_ringsD, bdy_ringsD, Dto0, points0] = quad_tessellation(width, height, plot_toggle);

% find pointsD which have approximated to boundary of target shape for instance, 
% pointsD with largest y-value can be mapped to upper part of target shape, sayig circle

boundB = 1:2*width;
boundL = 1+2*width: 4*width+1: (8*width+2)*height+2*width;
boundR = 6*width+1: 4*width+1: (8*width+2)*height+2*width;
boundT = boundL(end)+1:2:(8*width+2)*height+2*width;

%% find actual boundary and constuct initial guess in the deplyed space
ymin = min(points0(:,2));
xval = sort(unique(points0(:,1)));
xmax = xval(end);
boundRD = find(points0(Dto0,1) == xmax);
boundBD = find(points0(Dto0,2) == ymin);

% find the boundary edges in edgesD
edges_bottom = edgesD(ismember(edgesD(:,1),boundBD) & ismember(edgesD(:,2),boundBD),:);
edges_right = edgesD(ismember(edgesD(:,1),boundRD) & ismember(edgesD(:,2),boundRD),:);

% Load the target shape
shape = str2func(shape_name);  
[spline_boundR, spline_boundT, spline_boundL, spline_boundB] = shape();

% construct the initial guess in the deployed space
pointsD = compute_initial_map(pointsD_standard, shape_name, initial_map_type, ...
    scale_factor, boundR, boundT, boundL, boundB);

% plot the initial guess
figure(4)
clf
axis equal
axis off
hold on
plot_faces_generic(pointsD, face_setsD, 4)

plot(pointsD(boundR,1), pointsD(boundR,2), 'or')
plot(pointsD(boundT,1), pointsD(boundT,2), 'og')
plot(pointsD(boundL,1), pointsD(boundL,2), 'ob')
plot(pointsD(boundB,1), pointsD(boundB,2), 'oy')

fnplt(spline_boundR, [0 1], 'r', .5)
fnplt(spline_boundT, [0 1], 'g', .5)
fnplt(spline_boundL, [0 1], 'b', .5)
fnplt(spline_boundB, [0 1], 'y', .5)

title('Initial guess in the deployed space');
%% Constrained optimization

% optimization setup

% for the objective function
same_face_adjs = find_smooth_faces(width, height);
same_face_adjs = {same_face_adjs}; % for compatibility with OBJ

boundary_nodes_cell = {boundR, boundT, boundL, boundB};
boundary_target_splines_cell = {spline_boundR, spline_boundT, spline_boundL, spline_boundB};

options = optimoptions(@fmincon, ...
    'Display', 'iter-detailed', ...
    'Algorithm', 'interior-point', ...  % 'sqp' or 'interior-point'
    'SpecifyObjectiveGradient', true, ... % please keep it as true
    'SpecifyConstraintGradient', true, ... % please keep it as true
    'MaxFunctionEvaluations', 3000, ... % or try a larger number
    'MaxIterations', 200, ... % or try a larger number
    'ConstraintTolerance', 1e-6, ... % or try a smaller number
    'StepTolerance', 1e-6, ... % or try a smaller number
    'ScaleProblem', 'obj-and-constr', ... 
    'PlotFcn', {@optimplotfval, @optimplotconstrviolation,@optimplotfirstorderopt});

%% main optimization procedure
tic;
[solved_pointsD, ~, ~, ~] = fmincon(@(x)OBJ_regularization( ...
    decompose_v(x), face_setsD, same_face_adjs), ... objective function
    compose_v(pointsD), ... initial point
    [], [], [], [], [], [], ... linear constraints
    @(x)all_constraint_residual_and_jacobian( ... 
    decompose_v(x), ... initial point
    edgesD, edge_pairsD, anglesD, int_ringsD, bdy_ringsD, ... stencils
    boundary_nodes_cell, ... boundary nodes
    boundary_target_splines_cell, ... target shape
    [], [], ... non-overlap
    rectangular_ratio, edges_bottom, edges_right), ... contracted boundary shape control
    options ... optimization options
    );
toc;

% decompose the solved deployed structure
solved_pointsD = decompose_v(solved_pointsD);

%% plot the results

% plot the optimized deployed structure
h = figure(5);
clf
hold on
axis off
axis equal

bdy_ptD = [boundB, boundR, boundT, boundL];
plot(solved_pointsD(bdy_ptD, 1), solved_pointsD(bdy_ptD, 2), 'or');
plot_faces_generic(solved_pointsD, face_setsD, h)

before = findall(gca);
fnplt(spline_boundR, [0 1], 'k', 2)
added = setdiff(findall(gca), before);
set(added, 'Color', [201 0 22 200]/255);

before = findall(gca);
fnplt(spline_boundT, [0 1], 'k', 2)
added = setdiff(findall(gca), before);
set(added, 'Color', [201 0 22 200]/255);

before = findall(gca);
fnplt(spline_boundL, [0 1], 'k', 2)
added = setdiff(findall(gca), before);
set(added, 'Color', [201 0 22 200]/255);

before = findall(gca);
fnplt(spline_boundB, [0 1], 'k', 2)
added = setdiff(findall(gca), before);
set(added, 'Color', [201 0 22 200]/255);

title('Optimized pattern in the deployed space');


