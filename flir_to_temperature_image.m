function T = flir_to_temperature_image(filename)

setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%% Read metadata from JPG - exiftool application
% filename = 'example.jpg';
[~, img_data_json] = system(['exiftool -b -j ', filename]);
img_data = jsondecode(img_data_json);

%% Read RAW data - first 7 characters are encode information (base64:)
input = img_data.RawThermalImage(8:end);
%% Decode image from bytes 
if ischar(input) 
    input = uint8(input); 
end

output = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(input), 'uint16')';

output2 = output(103:end);
I1 = reshape(output2, [320 240]);
% I = I1''';
I = rot90(I1',1);

% imshow(I,[]);
%% Read information from metadata
tAtmC = str2double(img_data.AtmosphericTemperature(1:end-2));
Air_Temp = tAtmC + 273.15;
humidity = str2double(img_data.RelativeHumidity(1:end-2))/100;
Emissivity = img_data.Emissivity;
Refl_Temp = str2double(img_data.ReflectedApparentTemperature(1:end-2))+ 273.15;
Distance = str2double(img_data.ObjectDistance(1:end-2));

Alpha1 = img_data.AtmosphericTransAlpha1;
Alpha2 = img_data.AtmosphericTransAlpha2;
Beta1 = img_data.AtmosphericTransBeta1;
Beta2 = img_data.AtmosphericTransBeta2;
X = img_data.AtmosphericTransX;

S = double(I);
PlanckR1 = img_data.PlanckR1;
PlanckR2 = img_data.PlanckR2;
PlanckB = img_data.PlanckB;
PlanckO = img_data.PlanckO;
PlanckF = img_data.PlanckF;

%% Compute some parameters
h2o = humidity*exp(1.5587+0.06939*tAtmC-0.00027816*tAtmC^2+0.00000068455*tAtmC^3);
tau = X*exp(-sqrt(Distance)*(Alpha1+Beta1*sqrt(h2o)))+(1-X)*exp(-sqrt(Distance)*(Alpha2+Beta2*sqrt(h2o)));

Raw_Atm = PlanckR1/(PlanckR2*(exp(PlanckB/(Air_Temp))-PlanckF))-PlanckO;
Raw_Atm_tau = Raw_Atm*(1 - tau);

Raw_Refl = PlanckR1/(PlanckR2*(exp(PlanckB/(Refl_Temp))-PlanckF))-PlanckO;
Raw_Refl_em = Raw_Refl*(1-Emissivity)*tau;

Raw_Obj = (S-Raw_Refl_em-Raw_Atm_tau)/Emissivity/tau;

%% Convert RAW data to temperatures
T = PlanckB./log(PlanckR1./(PlanckR2*(Raw_Obj+PlanckO))+PlanckF)-273.15;
end

