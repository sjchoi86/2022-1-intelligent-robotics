%
% MATLAB tutorial
%
ccc
%% CCC
ccc % abbreviation for clc, clear all, and close all

%% Variables
ccc
a = 1;
b = 'b';
c = 'abc';
d = [1,2;3,4];

%% Print
ccc
fprintf(2,'Hello, World.\n');
fprintf("Hello, World.\n"); % bot works
disp('Hello, World.');
str = 'Hello, World.';
fprintf("%s [%d] [%.3f] \n",str,123,123.4567);

%% Simple plot
ccc
figure(1); hold on;
plot(rand(10,1),'o-','color','r','markersize',10,'linewidth',2);
plot(rand(10,1),'^-','color','b','markersize',10,'linewidth',2);
title('Plot Lines','fontsize',15);
xlabel('X'); ylabel('Y');
grid on; drawnow;

%% That's it. For the rest of functionalities, check out yourself.
ccc
