function [ null1,null2 ] = NullclinesFromXPP( filename )
%Loads nullclines from an XPP .dat file.
%%
nullclines = importdata(filename);

null1 = [nullclines.data(nullclines.data(:,3)==1,1),nullclines.data(nullclines.data(:,3)==1,2)];
null2 = [nullclines.data(nullclines.data(:,3)==2,1),nullclines.data(nullclines.data(:,3)==2,2)];

end

