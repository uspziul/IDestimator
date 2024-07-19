function [pxFiltrado,pyFiltrado] = IDestimator(matriz,plota,color,thickness)
%Matriz is a dataset in which the rows are the samples and the collumns are
%the dimensions. This function will create a vector log2(r) named px, and a
%vector log2(C(r)) named py, these axis are ploted and the slope in
%selected regions will be shown if plota==1 (by default plota=1),
%otherwise, vectors px and py will be returned to the workspace.

%This function is based on the method proposed by Grassberger and Procaccia
%(1983) to estimate intrinsic dimensionality (ID) and adapted by Montalvão
%et al. (2020) in which to get the ID, the number of samples C(r) that fall
%inside a hyper-cube of edge r is counted as r changes, the ID of the data
%is found by the slope of the curve on the plot log2(r) versus log2(C(r)).
%The edge size depends on the density of data and the scale of the data
%(usually found by the differential entropy of the data divided by the
%number of samples). Refer to Campadelli (2015) for more information on how
%to select r. This function also includes a bias compensation method
%developed in Montalvão et al. (2020).

if nargin<2
    plota=1;
end

% Flag to show or not the text representing the slop over the graph
showSlope = 1;

% Amount of text shown in the plot
nTexto = 20;


% Computes the intersample supremum norm
r = pdist(matriz,'chebychev');

r=r*2;
rSorted=sort(r);
lr = length(r);

py = log2((1:lr)/lr);
px = log2(rSorted);

% Remove duplicates and values to close to each other
[pxFiltrado, ipx, ~] = unique(round(px,3));
pxFiltrado = pxFiltrado';
pyFiltrado = py(ipx)';

% Frees up space in the memory
clear px
clear py

% Set text position to show slope in the graph
[pxtexto,inds] = subamostra(pxFiltrado,20);
pxtexto = pxtexto';
pytexto = pyFiltrado(inds)';

% Initializes slope (as) and y-intersect (h)
as = zeros(size(inds));
h = zeros(size(inds));

centros = zeros(length(as),2);
tira = [];

% Loop to calculate the average slope shown by text
for i = 1 : length(inds)
    if inds(i) > nTexto && ((inds(i)+nTexto) < length(pxFiltrado))

        % Computes the average slope between points represented by the text
        % in the plot
        temp = pinv([pxFiltrado(inds(i)-nTexto:inds(i)+nTexto) ones(2*nTexto + 1,1)]) * pyFiltrado(inds(i)-nTexto:inds(i)+nTexto);
        
        % Stores the average slope in as
        as(i) = temp(1);

        % Stores the averare y-intersect in h
        h(i) = temp(2);

        centros(i,1) = pxtexto(i);
        centros(i,2) = pytexto(i);

    else
        tira = [tira; i];
    end
end
centros(tira,:) = [];
as(tira) = [];

% Plots the graph
if plota==1
    if nargin==4
        plot(pxFiltrado,pyFiltrado,color,'LineWidth',thickness);
    else
        plot(pxFiltrado,pyFiltrado);
    end
    if showSlope==1
        text(centros(:,1),centros(:,2),cellstr(num2str(as')));
        idmedcentro = mean(as(round(0.333*length(as):round(0.667*length(as)))));
    end
    title(idmedcentro);
    grid on 
end

% ----------------------Bias correction part-------------------------------

% Implementation of the bias compensation method developed in Montalvão et
% al. (2020)

d = idmedcentro; %prototype;
dTeste = d:.1:d+8;
dTeste(dTeste==0)=[];

[li,~]=size(matriz);

N = li;
rTraco = 1./(1+N.^(1./dTeste));
d0 = dTeste.*(1-(rTraco./(2-rTraco)));
[~,pos] = min(abs(dTeste-d0));
id = dTeste(pos);

for k=1:length(dTeste)
    rb=1/(1+N^(1/dTeste(k)));
    d0(k)=dTeste(k)*(1-rb/(2-rb));
    dh(k)=dTeste(k)*(rb/(2-rb)*log2(rb)+log2(2-rb));
end

[~,menor] = min(abs(d0 - d));
idmedcentro = dTeste(menor);
title(strcat("ID before bias compensation = ", num2str(d),". ID after bias compensation = ", num2str(idmedcentro)));
