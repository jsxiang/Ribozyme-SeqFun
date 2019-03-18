%% Demo for Standard function.
%
%   This script performs demo of the Standard function.
%   
%   Database is loaded and two sets are crated. The first set is
%   standardized and its Mean and Std are used on second set to be
%   standardized.
%

%% Load the database.
%
%   Loading the needed database for the purpose of this demo. In this demo
%   the built-in database acetylene.mat is used.
%

load acetylene.mat

%% Prepare the sets.
%
%   Prepare the data sets. Divide the data in two sets.
%

x = [ x1, x2, x3, y ];         %   Prepare the data.

Set1    = x ( 1 : 8, : );     %   Set 1 

Set2    = x ( 9 : end, : );   %   Set 2

%% Standardize Set 1
%
%   Use Standard function on Set 1
%   Standardize by Row.
%

[ SetStandard1, CMean, CStd] = Standard(Set1, 1);

%% Standardize Set 2
%
%   Use Standard function on Set2. 
%   Use input CMean and CStd
%   Standardize by Row.

[ SetStandard2] = Standard(Set2, 1,CMean,CStd);
