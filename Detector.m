classdef Detector < handle
 %Detector ������������ ����������� �����

 properties
 %���������� ��������� ���������
 thr %����� �����������
 Fa %����������� ������ ������

 detections %�������� ������ ������
 end

 methods
 function obj = Detector(fa)
 obj.Fa = fa;
 end
 function detect(obj,inData,r,a) %����������� �����
 %��������� ����� � �������� ������� ��������� �� ����� r - a
 temp = abs(inData(r,a));
 s2 = var(temp(:));
 obj.thr = 2.5*sqrt(-2*log(obj.Fa)*s2);
 obj.detections = inData > obj.thr;
 end
 function show(obj,hInData)
 figure;
 %����� ������
 subplot(1,2,1);
 temp = abs(hInData(:));
 q = quantile(temp,0.9);
 imagesc(abs(hInData),[0 q]);
 xlabel('����, ���');
 ylabel('���������, ���');

 if ~isempty(obj.detections)
 subplot(1,2,2);
 %������ ����� ������
 imagesc(abs(obj.detections),[0 1]);
 xlabel('����, ���');
 ylabel('���������, ���');
 end
 end
 end
end