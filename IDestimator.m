function [pxFiltrado,pyFiltrado] = IDestimator(matriz,plota)
%Matriz is a dataset in which the rows are the samples and the collumns are
%the dimensions. This function will create a vector log2(r) named px, and a
%vector log2(C(r)) named py, these axis are ploted and the slope in
%selected regions will be shown if plota==1 (by default plota=1), 
%otherwise, vectors px and py will be returned to the workspace. 

%This function is based on the method proposed by Grassberger and
%Procaccia (1983) to estimate ID and adapted by Montalvão et al. (2020) in
%which To get the ID, the number of samples C(r) that fall inside a 
%hyper-cube of edge r is counted as r changes, the ID of the data is found 
%by the slope of the curve on the plot log2(r) versus log2(C(r)). To get
%the ID you have to choose an edge size r returned by 2^px and find the
%slope of the plot at this value. The edge size depends on the density of
%data and the scale of the data (usually found by the differential entropy
%of the data divided by the number of samples). Refer to Campadelli (2015)
%for more information on how to select r.

if nargin<2
    plota=1;
end

[li,co] = size(matriz);

x=matriz;

r = pdist(matriz,'chebychev');

r=r*2;
rSorted=sort(r);
lr = length(r);

py = log2((1:lr)/lr);
px = log2(rSorted);

pxFiltrado=px(1);
pyFiltrado=py(1);
for i = 2:length(px)
    if px(i)~=px(i-1)
        pxFiltrado=[pxFiltrado;px(i)];
        pyFiltrado=[pyFiltrado;py(i)];
    end
end
clear px
clear py

[pxtexto,inds] = subamostra(pxFiltrado,20);
pxtexto = pxtexto';
pytexto = pyFiltrado(inds);
dydx = diff(pytexto)./diff(pxtexto);

if plota==1
    plot(pxFiltrado,pyFiltrado);
    axis equal
    text(pxtexto(1:end-1),pytexto(1:end-1),cellstr(num2str(dydx)));
    grid on 
end
