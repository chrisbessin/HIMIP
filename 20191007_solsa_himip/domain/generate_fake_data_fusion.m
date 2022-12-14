function generate_fake_data_fusion(size_XY,size_VNIR,size_SWIR,dir_data)
%
%
size_X = size_XY(1);
size_Y = size_XY(2);
size_VNIR_X = size_VNIR(1);
size_VNIR_Y = size_VNIR(2);
size_SWIR_X = size_SWIR(1);
size_SWIR_Y = size_SWIR(2);

% XYZ
def_x = linspace(-25,25,size_X);
X = repmat(def_x',1,size_Y);
def_y = linspace(-5,50,size_Y);
Y = repmat(def_y,size_X,1);
Z = -0.25*X.^2 + 160;

% RGB
def_x = linspace(1,size_X,size_X);
index_x_rgb = repmat(def_x',1,size_Y);
def_y = linspace(1,size_Y,size_Y);
index_y_rgb = repmat(def_y,size_X,1);

% VNIR
def_x = linspace(1,size_VNIR_X,size_X);
index_x_vnir = repmat(round(def_x',0),1,size_Y);
def_y = linspace(1,size_VNIR_Y,size_Y);
index_y_vnir = repmat(round(def_y,0),size_X,1);

% SWIR
def_x = linspace(1,size_SWIR_X,size_X);
index_x_swir = repmat(round(def_x',0),1,size_Y);
def_y = linspace(1,size_SWIR_Y,size_Y);
index_y_swir = repmat(round(def_y,0),size_X,1);

ref_dir = pwd;
cd(dir_data)
save('XYZ','X','Y','Z')

cd('RGB')
save('index_x_rgb+index_y_rgb','index_x_rgb','index_y_rgb')

cd('..\VNIR')
save('index_x_vnir+index_y_vnir','index_x_vnir','index_y_vnir')

cd('..\SWIR')
save('index_x_swir+index_y_swir','index_x_swir','index_y_swir')


cd(ref_dir)
end

