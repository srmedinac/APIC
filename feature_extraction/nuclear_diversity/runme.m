
clear all;

%% step 1: load image
L = imread("C:\Users\smedin7\OneDrive - Emory University\Documents\data\renal_project_liping\nuclei_segmentation\S20-18383-A6 - 2022-06-02 12.28.56\x6144_y57344_w2048_h2048.png");
img = imread("C:\Users\smedin7\OneDrive - Emory University\Documents\data\renal_project_liping\patches\S20-18383-A6 - 2022-06-02 12.28.56\tiles\x6144_y57344_w2048_h2048.jpeg");

%% step 2: segment nuclei and save boundaries
[properties, bounds] = LNuclear_regionProperties(img,L);

ctemp=[properties.Centroid];
bounds.centroid_c=ctemp(1:2:end);
bounds.centroid_r=ctemp(2:2:end);
%% step 3: extract features
% cellular diversity features (1*20*13*6=1560)
para.CGalpha_min=0.44; para.CGalpha_max=0.44;% larger the alpha, sparse the local cell graph
para.alpha_res=0.02;
para.radius=0.2;

CGinfo=[]; % cell graph information
feats=[];
feats_description=[];
set_alpha=[para.CGalpha_min:para.alpha_res:para.CGalpha_max];
for f=1:length(set_alpha)
    curpara.alpha=set_alpha(f); curpara.radius=para.radius;
    [feat,feat_description,CGinfo{f}] = Lextract_CellularDiversity_features(bounds,properties,...
        curpara.alpha,curpara.radius);
    %%% get feature description
    temp= feat_description;
    for i=1:length(temp)
        cur=temp{i};
        str=sprintf('-a=%.2f',curpara.alpha);
        cur(end+1:end+length(str))=str;
        temp{i}=cur;
    end

    feats=cat(2,feats,feat); 
    feats_description=cat(2,feats_description,temp);
end
disp('feature extraction done. please check the variable feats and feats_description for cellular diversity features and their feature names');