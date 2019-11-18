classdef ConvData < handle
%ConvData ������������ ������������� ���������� ������

properties
m_RawData %��������� �� ������ ������ RawData

%���������� ��������� �������� ��� �������
ski_tau = 20e-9; %������������ ��������, ���
ski_amp = 1; %��������� �������� ��������, ���

%���������� ��������� �������� ��� �������
lfm_tau = 5e-6; %������������ ��������, ���
lfm_df = 50e6; %�������� ������� �������, ��
lfm_amp = 1; %��������� �������� ��������, ���

%���������� ��������� ������������� �������
f_samp = 400e6; %��
f_if = 0; %������������� ������� ������� ����� �����������������

h_ski %������� ��� ������
h_lfm %������� ��� ������
convolved_data %�������� ������ ������

%��������� ����������� �������
med %������� �������� ������ �� �����
end

methods
function obj = ConvData(hRawData) %����������� �������
obj.m_RawData = hRawData; %������� ����� � �������� ������ RawData
%������� ������� �������
time = 0:1/obj.f_samp:obj.ski_tau-1/obj.f_samp; %������� ������� �� ����� ���������
%������� ������ �� ������� (5)
obj.h_ski = obj.ski_amp*exp(1j*2*pi*obj.f_if*time);
obj.h_ski = fliplr(obj.h_ski); %� ����. ������� �� ��� ���������� ��������� ������� �� �������

time = 0:1/obj.f_samp:obj.lfm_tau-1/obj.f_samp; %������� ������� �� ����� ���������
%������� ������ �� ������� (6)
obj.h_lfm = obj.lfm_amp*exp(1j*2*pi*(obj.f_if+obj.lfm_df/obj.lfm_tau*time).*time);
obj.h_lfm = fliplr(obj.h_lfm); %� ����. ������� �� ��� ���������� ��������� ������� �� �������

%���������� ��������� �����������
obj.med = median(abs(obj.m_RawData.raw(:)));
end
function convole_time(obj,type) %������� �� ��������� �������
obj.convolved_data = []; %������� ���������� ������� ����������
temp = obj.m_RawData.raw; %�������� ����� ��� �������
n = size(temp,2); %����� ������� �� ���������
if strcmp(type,'ski') %������������ ���
m = length(obj.h_ski);
obj.convolved_data = zeros(m+size(temp,1)-1,n);
for ii = 1:n
obj.convolved_data(:,ii) = conv(temp(:,ii),obj.h_ski);
end
else %������������ ���
m = length(obj.h_lfm);
obj.convolved_data = zeros(m+size(temp,1)-1,n);
for ii = 1:n
obj.convolved_data(:,ii) = conv(temp(:,ii),obj.h_lfm);
end
end
obj.show;
end
function convole_fft(obj,type) %������� � ��������� �������
obj.convolved_data = []; %������� ���������� ������� ����������
temp = obj.m_RawData.raw; %�������� ����� ��� �������
n = size(temp,2); %����� ������� �� ���������
if strcmp(type,'ski') %������������ ���
m = length(obj.h_ski);
%������� ������� ���� �� ������� (10)
wnd = hamming(m)';
%������� ������� ������ ������� (�.�. �� �������� ��
%�������, �� ����������� ���������� �� �����), �����
%������� ������ ���� ������ �������� ������� ������� �
%������� ������� (������� 7)
H_ski = fft(wnd.*obj.h_ski,m+size(temp,1))';
%���������� ������ �� ���� ������������ �������
H_ski = repmat(H_ski,1,n);
 
%��������� ������ ������� � ������ ������� ������ (�.8)
H_RawData = fft(temp,m+size(temp,1),1);
%������������ ������� ������� (�.9)
obj.convolved_data = ifft(H_RawData.*H_ski);
else %������������ ���
m = length(obj.h_lfm);
%������� ������� ���� �� ������� (10)
wnd = hamming(m)';
%������� ������� ������ ������� (�.�. �� �������� ��
%�������, �� ����������� ���������� �� �����), �����
%������� ������ ���� ������ �������� ������� ������� �
%������� ������� �.(7)
H_lfm = fft(wnd.*obj.h_lfm,m+size(temp,1))';
%���������� ������ �� ���� ������������ �������
H_lfm = repmat(H_lfm,1,n);
%��������� ������ ������� � ������ ������� ������ �.(8)
H_RawData = fft(temp,m+size(temp,1),1);
%������������ ������� ������� �.(9)
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