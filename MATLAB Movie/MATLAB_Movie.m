% This is for .xyz coordinates file format with only sphere for representing atoms
clc
clear all
%% Rough
clc
s = 1;
rng(s);
IE = rand(10,1)*100;
dummy_IE = zeros(10,1);
Order = zeros(10,1);
dummy_IE(1) = IE(1);
for i = 1:10
    dummy_IE(i) = IE(i);
    Order(i) = i;
    for j = 1:i
        if IE(i) > dummy_IE(j)
            for m = i:-1:j+1
                dummy_IE(m) = dummy_IE(m-1);
                Order(m) = Order(m-1);
            end
            dummy_IE(j) = IE(i);
            Order(j) = i;
            break;
        end
    end
    [IE dummy_IE]  
end
% [IE dummy_IE]    
%% Read lines of the text file
A = readlines("Coordinates_short0.xyz");
NFrames = 0;
NMolec = [];
NMetal = 144;
CoordMol = [];
Molec = [];
Frame = [];
Metals = ["Fe", "Pt", "Ru", "Cu"];
MetalColor = [0.8500 0.3250 0.0980; 0.8157 0.8157 0.8784; 0.2 0.5 0; 0.7843 0.502 0.2];
MetalCoord = [];
MetalChange = [];
bool_MetalChange = 0;
% Extract all the data and store in arrays
i = 1;
nmet = 0;
nmol = 0;
NFrames = 0;
First = 1;
read = 0;
MetalType = 0;
while i<length(A) && A(i) ~= "END"
    disp("Done reading line: " + i);
    dum2 = A(i); % Number of Atoms in the frame
    %Extracting the first word from the line:
    NAtoms = str2num(dum2); %Number of atoms in the frame
    NAtoms = NAtoms - NMetal;
    i = i + 1;
    dum2 = A(i);
    NFrames = NFrames + 1;
    if NFrames > 0
        read = 1;
        Frame = [Frame; str2num(dum2)]; %Frame Number
        NMolec = [NMolec; NAtoms];
    end
    i = i + 1;
    if First == 1 %Collect the Metal atom coordinate information
        for Metal_line = 1:NMetal
            dum2 = A(i);
            if read == 1
                MetalCoord = [MetalCoord; ReadAtoms(dum2)];
            end
            if MetalType == 0
                for met = 1:length(Metals)
                    if Metals(met) == MetalCoord(1,1)
                        MetalType = met;
                    end
                end
            end
            i = i + 1;
        end
        First = 0;
    else
        i = i + NMetal;
    end
    for Molec_line = 1:NAtoms
        dum2 = A(i);
        if read == 1
            CoordMol = [CoordMol; ReadAtoms(dum2)];
        end
        i = i + 1;
    end
end
%% Plotting metal and adsorbate atoms
SPHERE_RES = 20; % resolution for your spheres
Metal_Sphere = [1.26, 1.387, 1.387, 1.387]; %Radius of the spheres of Fe, Pt, Ru, Cu
[xb, yb, zb] = sphere(SPHERE_RES);
axis equal;
set(gcf,'Color','w');
metal_plot = plot(1,1); %Initializing a plot variable
init3dview();
% Plotting Metal atoms which remain on the figure throughout
% delete(metal_plot);
[Metal, MCoord] = plotAtoms(NMetal, MetalCoord);
material([0 1 0]);
v = MCoord;
for i=1:NMetal
    for k=1:length(Metals)
        if Metal(i) == Metals(k)
            SPHERE_RAD = Metal_Sphere(k);
            col = MetalColor(MetalType, :);
            metal_plot = [metal_plot; surface(SPHERE_RAD*xb+v(i,1), SPHERE_RAD*yb+v(i,2), SPHERE_RAD*zb+v(i,3), 'facecolor', col, 'EdgeColor','none')]; % This line is changed
        end
    end
end
% Plotting Adsorbates
molec_plot = [];
SPHERE = [0.8, 0.53, 0.8, 0.8]; % radius of spheres in the order N H C O
M(NFrames) = struct('cdata',[],'colormap',[]);
MolecPlotted = 0;
for j=1:NFrames
    delete(molec_plot);
    molec_plot = [];
    if j==1
        CF = CoordMol(1:NMolec(j),:);
    else
        CF = CoordMol(MolecPlotted+1:MolecPlotted + NMolec(j),:);
    end
    [Molec, v] = plotAtoms(NMolec(j), CF);
    MolecPlotted = MolecPlotted + NMolec(j);
    for i = 1:NMolec(j)
        if Molec(i) == "H"
            col = [1, 1, 1];
            SPHERE_RAD = SPHERE(2);
            molec_plot = [molec_plot; surface(SPHERE_RAD*xb+v(i,1), SPHERE_RAD*yb+v(i,2), SPHERE_RAD*zb+v(i,3), 'facecolor', col, 'EdgeColor','none')]; % This line is changed
        elseif Molec(i) == "N"
            col = "b";
            SPHERE_RAD = SPHERE(1);
            molec_plot = [molec_plot; surface(SPHERE_RAD*xb+v(i,1), SPHERE_RAD*yb+v(i,2), SPHERE_RAD*zb+v(i,3), 'facecolor', col, 'EdgeColor','none')]; % This line is changed
        elseif Molec(i) == "C"
            col = [1, 1, 1]*0.5647;
            SPHERE_RAD = SPHERE(3);
            molec_plot = [molec_plot; surface(SPHERE_RAD*xb+v(i,1), SPHERE_RAD*yb+v(i,2), SPHERE_RAD*zb+v(i,3), 'facecolor', col, 'EdgeColor','none')]; % This line is changed
        elseif Molec(i) == "O"
            col = [1, 0.051, 0.051];
            SPHERE_RAD = SPHERE(4);
            molec_plot = [molec_plot; surface(SPHERE_RAD*xb+v(i,1), SPHERE_RAD*yb+v(i,2), SPHERE_RAD*zb+v(i,3), 'facecolor', col, 'EdgeColor','none')]; % This line is changed
        end
    end
    M(j) = getframe;
end
%%
v = VideoWriter('movie.avi');
v.FrameRate = 2;
open(v);
for k = 1:NFrames 
   writeVideo(v,M(k));
end
close(v);
%% Rough for checking
% clc
% word = Line2word(A(1));
%% Functions
% Get all the words in an array
function words = Line2word(line)
    words = [];
    j = 1;
    N_line = strlength(line);
    while j<=N_line+1
        dum1 = "";
        while extract(line,j)~=" "
            dum1 = dum1 + extract(line,j);
            j = j + 1;
            if j>N_line
                break;
            end
        end
        words = [words dum1];
        if j>N_line
                break;
        end
        while extract(line,j)==" "
            j = j + 1;
            if j>N_line
                break;
            end
        end  
        if j>N_line
                break;
        end
    end
end
function Coords = ReadAtoms(A)
    Coords = [];
    NumbChar = strlength(A);
    i = 1;
    while i<=NumbChar
        dum1 = "";
        word = "";
        while dum1 ~= " "
            if (i<=NumbChar)
                word = word + dum1;
                dum1 = extract(A,i);
            else
                word = word + dum1;
                break;
            end
            i = i + 1;
        end
        Coords = [Coords string(word)];
    end
end
function boo = Contained(str, arra)
    boo = false;
    for i =1:length(arra)
%         arra(i)
        if str==arra(i)
            boo = true;
        end
    end
end
function [Metal, MCoord] = plotAtoms(NMet, Coords)
    Metal = [];
    MCoord = [];
    for i = 1:NMet
        A = Coords(i,:);
        Metal = [Metal; A(1)];
        MCoord = [MCoord; str2double(A(2)), str2double(A(3)), str2double(A(4))];
    end
end
function init3dview()
  axis vis3d
  axis off;
  xlim([-5 45])
  ylim([-10 45])
  zlim([0 10])
%   view(3)
  camlight
end