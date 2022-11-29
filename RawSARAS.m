%% SARAS (Simulador Académico de Redes Acústicas Subacuáticas)
% Autor: Angel Velasco      % Última actualización: 29/11/2022
clc
clear

%% Variables
dataset.seed =6;
dataset.nodeNo = 25;
dataset.area = [500 250 0]; % [m]
dataset.range = 100;        % [m]
dataset.nodeSrc = 1;
dataset.nodeDst = 25;
dataset.frecuency = 10;     % [KHz]
dataset.power = 2;          % [W]
dataset.data = 5000;        % [b]

%% RNG
if dataset.seed ~= 0
    rng(dataset.seed);
end

%% Crear tabla de nodos
varNames = {'X','Y','Naturaleza'};
varTypes = {'single','single','categorical'};
dataset.nodePosition = table('Size',[dataset.nodeNo 3], ...
    'VariableTypes',varTypes,'VariableNames',varNames);
dataset.nodePosition.Properties.VariableUnits = {'m','m',''};
dataset.nodePosition.Naturaleza = categorical( ...
    dataset.nodePosition.Naturaleza,{'Ordinario','Sumidero'});
dataset.nodePosition.Naturaleza = fillmissing( ...
    dataset.nodePosition.Naturaleza, 'constant', "Ordinario");

%% Generación de nodos
for i = 1 : dataset.nodeNo
    dataset.nodePosition.X(i) = randi([1 dataset.area(1)]);
    dataset.nodePosition.Y(i) = randi([1 dataset.area(2)]);
end
disp(dataset.nodePosition)

%% Cálculo de distancia euclidiana
for i = 1 : dataset.nodeNo
    for j = 1: dataset.nodeNo
        garbage.x1 = dataset.nodePosition.X(i);
        garbage.x2 = dataset.nodePosition.X(j);
        garbage.y1 = dataset.nodePosition.Y(i);
        garbage.y2 = dataset.nodePosition.Y(j);
        dataset.euclidiana(i,j) = sqrt((garbage.x1 - ...
            garbage.x2) ^2 + (garbage.y1 - garbage.y2)^2);
    end
end

%% Construcción de grafo
dataset.weights = lt(dataset.euclidiana,dataset.range);
G=graph(dataset.weights,'omitselfloops');
for a = 1 : height(G.Edges)
    garbage.s = G.Edges.EndNodes(a,1);
    garbage.t = G.Edges.EndNodes(a,2);
    garbage.Z(a,:) = dataset.euclidiana(garbage.s,garbage.t);
end
G.Edges.Euclidiana = garbage.Z(:,1);

%% Ploteo
close all
figure('units','normalized','innerposition',[0 0 1 1], ...
    'MenuBar','none')
p = plot(G,'XData',(dataset.nodePosition.X),'YData',( ...
    dataset.nodePosition.Y),'MarkerSize',8,'NodeColor', ...
    '0.80,0.80,0.80','EdgeColor','0.65,0.65,0.65', ...
    'LineWidth',0.7,'NodeLabelColor','0.90,0.90,0.90');
line(dataset.nodePosition.X(dataset.nodeSrc), ...
    dataset.nodePosition.Y(dataset.nodeSrc),'color','green', ...
    'marker','o','linestyle','none','markersize',50)
line(dataset.nodePosition.X(dataset.nodeDst), ...
    dataset.nodePosition.Y(dataset.nodeDst),'color','green', ...
    'marker','o','linestyle','none','markersize',50)
grid on
hold on
title('Dispersión de Nodos')

%% Enrutamiento
Gpath = shortestpathtree(G,dataset.nodeSrc,dataset.nodeDst);
dataset.routepath = shortestpath(G,dataset.nodeSrc, ...
    dataset.nodeDst);
dataset.hopsNo = length(dataset.routepath)-1;
highlight(p,dataset.routepath,'NodeColor','green', ...
    'NodeLabelColor','green')
disp(dataset.routepath)
disp(dataset.hopsNo)

%% Transmisión de datos
f = dataset.frecuency;
p = dataset.power;
SNRsum = 0;
for i = 1:dataset.hopsNo
    d = dataset.euclidiana(dataset.routepath(i), ...
        dataset.routepath(i+1));
    absortC = (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) ...
        + 2.75e-4*f^2 + 0.003; % Coeficiente de absorción
    SS = 20*log10(d); % Propagación esférica
    Slevel = 10*(log10(p) - log10(4*pi*d^2) ...
        - log10(0.67e-18));
    Tloss = SS + absortC*d*10^(-3);
    Nlevel = 50 - 18*log10(f);
    Dindex = 0; % Debido a la transmisión omnidireccional
    aux = 10^(Slevel/10) - 10^(Tloss/10) - 10^(Nlevel/10) ...
        + 10^(Dindex/10);
    SNRsum = SNRsum + aux;
end
dataset.SNR = 10 * log10(SNRsum/dataset.hopsNo);
dataset.BER = 0.5 * (1 - sqrt((10^(dataset.SNR/10))/(1 + ...
    10^(dataset.SNR/10))));
dataset.SDR = (1 - dataset.BER)^(dataset.data);
disp(dataset.SNR)
disp(dataset.BER)
disp(dataset.SDR)

