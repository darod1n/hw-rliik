classdef AKP < handle
 %AKP ������������ ��������������� �������� �����

 properties
 %���������� ��������� ���������������
 wnd = 30 %����� �������� ��� ������ �����. �����������
 g = 10 %����������� �����������

 th_shift = 5000; %�������� ��� ������� �������� 5000 = -180 ��
 k_back = 10^(-13.3/20); %����������� �������� ��������������� �������

 akp_data %�������� ������ ������
 akp_ang %�������� ������ ����� ����� ������
 end

 methods
 function compensate(obj,hInData,hInAng) %����������� ������
 obj.akp_data = []; %������� ���������� ������� ����������
 obj.akp_ang = [];

 %��������� ������ ����������� �� �����
 a = (hInAng - obj.th_shift)>1;
 [~,col] = find(a > 0,1,'first');
 for ii = col:length(hInAng)
 u0 = hInData(:,ii); %������� 20
 angIdx = ii-obj.th_shift:ii-1;
 angIdx(angIdx < 1) = [];
 u1 = sum(hInData(:,angIdx),2)*obj.k_back; %������� 21 � 24
 K = obj.g*sum(u0(end-obj.wnd:end).*u1(end-obj.wnd:end),1)/...
 (1+obj.g*sum(u1(end-obj.wnd:end).^2,1));%������� 23
 obj.akp_data = [obj.akp_data (u0 - K*u1)];%������� 22
 obj.akp_ang = [obj.akp_ang hInAng(ii)];
 end
 end

 function show(obj,hInData,hInAng)
 figure;
 %����� ������
 subplot(1,2,1);
 temp = abs(hInData(:));
 q = quantile(temp,0.9);
 idxAng = (hInAng - obj.th_shift) > 1;
 imagesc(abs(hInData(:,idxAng)),[0 q]);
 xlabel('����, ���');
 ylabel('���������, ���');

 subplot(1,2,2);
 %������ ����� �����������
 imagesc(abs(obj.akp_data),[0 q]);
 xlabel('����, ���');
 ylabel('���������, ���');
 end

 function Kp = calc_kp(obj,r,a,hInData)
 %������ ������������ ���������� � ��������� �����
 Kp = 20*log10(abs(obj.akp_out)/abs(hInData(r,a)));
 end
 end
end