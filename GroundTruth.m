clear all, close all, clc
%% Import the data

addpath(genpath('./Data'))
load('Indian_pines_gt.mat');

%% Plot the ground truth

plot_classes(indian_pines_gt)
title('Ground Truth','FontSize',14)

saveas(gcf,'./Images/Indian_Pines_GroundTruth.png')