classdef ConvData < handle
%ConvData осуществляет согласованную фильтрацию данных

properties
m_RawData %указатель на объект класса RawData

%определяем параметры опорного ски сигнала
ski_tau = 20e-9; %длительность импульса, сек
ski_amp = 1; %амплитуда опорного импульса, отс

%определяем параметры опорного лчм сигнала
lfm_tau = 5e-6; %длительность импульса, сек
lfm_df = 50e6; %девиация частоты сигнала, Гц
lfm_amp = 1; %амплитуда опорного импульса, отс

%определяем параметры дискретизации сигнала
f_samp = 400e6; %Гц
f_if = 0; %промежуточная частота сигнала после гетеродинирования

h_ski %опорный ски сигнал
h_lfm %опорный лчм сигнал
convolved_data %выходной массив данных

%параметры отображения сигнала
med %медиана амплитуд данных на входе
end

methods
function obj = ConvData(hRawData) %конструктор объекта
obj.m_RawData = hRawData; %создаем связь с объектом класса RawData
%создаем опорные сигналы
time = 0:1/obj.f_samp:obj.ski_tau-1/obj.f_samp; %отсчеты времени на длине импульсов
%опорный сигнал по формуле (5)
obj.h_ski = obj.ski_amp*exp(1j*2*pi*obj.f_if*time);
obj.h_ski = fliplr(obj.h_ski); %в согл. фильтре ИХ это зеркальное отражение сигнала во времени

time = 0:1/obj.f_samp:obj.lfm_tau-1/obj.f_samp; %отсчеты времени на длине импульсов
%опорный сигнал по формуле (6)
obj.h_lfm = obj.lfm_amp*exp(1j*2*pi*(obj.f_if+obj.lfm_df/obj.lfm_tau*time).*time);
obj.h_lfm = fliplr(obj.h_lfm); %в согл. фильтре ИХ это зеркальное отражение сигнала во времени

%определяем параметры отображения
obj.med = median(abs(obj.m_RawData.raw(:)));
end
function convole_time(obj,type) %свертка во временной области
obj.convolved_data = []; %удаляем результаты прошлых вычислений
temp = obj.m_RawData.raw; %исходные данны для свертки
n = size(temp,2); %число сверток по дальности
if strcmp(type,'ski') %обрабатываем ски
m = length(obj.h_ski);
obj.convolved_data = zeros(m+size(temp,1)-1,n);
for ii = 1:n
obj.convolved_data(:,ii) = conv(temp(:,ii),obj.h_ski);
end
else %обрабатываем лчм
m = length(obj.h_lfm);
obj.convolved_data = zeros(m+size(temp,1)-1,n);
for ii = 1:n
obj.convolved_data(:,ii) = conv(temp(:,ii),obj.h_lfm);
end
end
obj.show;
end
function convole_fft(obj,type) %свертка в частотной области
obj.convolved_data = []; %удаляем результаты прошлых вычислений
temp = obj.m_RawData.raw; %исходные данны для свертки
n = size(temp,2); %число сверток по дальности
if strcmp(type,'ski') %обрабатываем ски
m = length(obj.h_ski);
%создаем весовое окно по формуле (10)
wnd = hamming(m)';
%создаем опорный спектр сигнала (т.к. он повернут во
%времени, то комплесного сопряжения не нужно), длина
%спектра должна быть суммой размеров опорной выборки и
%спектра сигнала (формула 7)
H_ski = fft(wnd.*obj.h_ski,m+size(temp,1))';
%размножаем спектр по всем азимутальным каналам
H_ski = repmat(H_ski,1,n);
 
%вычисляем спектр сигнала в каждом угловом канале (ф.8)
H_RawData = fft(temp,m+size(temp,1),1);
%осуществляем свертку сигнала (ф.9)
obj.convolved_data = ifft(H_RawData.*H_ski);
else %обрабатываем лчм
m = length(obj.h_lfm);
%создаем весовое окно по формуле (10)
wnd = hamming(m)';
%создаем опорный спектр сигнала (т.к. он повернут во
%времени, то комплесного сопряжения не нужно), длина
%спектра должна быть суммой размеров опорной выборки и
%спектра сигнала ф.(7)
H_lfm = fft(wnd.*obj.h_lfm,m+size(temp,1))';
%размножаем спектр по всем азимутальным каналам
H_lfm = repmat(H_lfm,1,n);
%вычисляем спектр сигнала в каждом угловом канале ф.(8)
H_RawData = fft(temp,m+size(temp,1),1);
%осуществляем свертку сигнала ф.(9)
obj.convolved_data = ifft(H_RawData.*H_lfm);
end
obj.show;
end
function show(obj)
figure;
subplot(1,2,1);
imagesc(abs(obj.m_RawData.raw),[0 100*abs(obj.med)]);
subplot(1,2,2);
imagesc(abs(obj.convolved_data),[0 100*abs(obj.med)]);
end
end

end