classdef RawDataLfm < handle           %��������� ����������� ������ hanlde � �������
                                    %RawData ��������� �������� ������ ��� � �� "�������������"
    properties
        fname %��� �����
        fpath %���� � �����
        raw %������ �������� � ����������� ���������-����
        angles %������ ����� ����� � ��������
        m %������� ���������� �������� � �����
        flag_raf%���� �������������
        mean_raw = zeros(1,3); %��� �������� �����
        disp_re = zeros(1,3); %��������� ������������ ����� �����
        disp_im = zeros(1,3); %��������� ������ ����� �����

        %�������� ���� ��� ���������
        zone1 = [1:1000];
        zone2 = [1001:2000];
        zone3 = [2001:16000];
    end

    methods
        function obj = RawDataLfm(fpath,fn) %����������� ������
            if nargin < 1
                obj.open;
            else
                obj.open(fpath,fn);
            end
            %����������� ������ ��� ����, ���� ���������, ��� � �������
            %������� ���������� ������������ ������ ����
            obj.show;
        end
        function open(obj,fpath,fn) %�������� ����� ��� ���������
            if nargin < 1
            %����� ����� ��� ���������
 [obj.fname, obj.fpath] = uigetfile('*.mat');
 else
 obj.fpath = fpath;
 obj.fname = fn;
 end
 %���������� ����� � float
 temp = (load([obj.fpath obj.fname],'lfmraw'));
 obj.raw = single(temp.lfmraw);
 temp = (load([obj.fpath obj.fname],'lfmangles'));
 obj.angles = single(temp.lfmangles);
 %���������� ������� ��� ����������� ������
 obj.m = median(abs(obj.raw(:)));
 %����� ����� �������������
 obj.flag_raf = 0; %��������� �� �����������
 end
 function show(obj) %����������� ������
 figure;
 imagesc(abs(obj.raw),[0 2*obj.m]);
 end
 function rafinate(obj,m)
 %������� ������: ��� ������ ��������� � m ��������� �������
 %���������
 %% ���� 1
 temp = obj.raw(obj.zone1(end)-m:obj.zone1(end),:);
 %������� ������ ������� �������� (��������� ��� ���������)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %���������� ���. �������� �� ������� (1)
 obj.mean_raw(1) = mean(temp);
 %���������� ��������� �� ������� (2)-(3)
 obj.disp_re(1) = var(real(temp)); %��������� ������������ ����� �����
 obj.disp_im(1) = var(imag(temp)); %��������� ������ ����� �����
 %����������� ��������� �� ������� (4)
 obj.raw(obj.zone1,:) = 1/sqrt(obj.disp_re(1))*(real(obj.raw(obj.zone1,:))-real(obj.mean_raw(1)))
+...
 1j/sqrt(obj.disp_im(1))*(imag(obj.raw(obj.zone1,:))-imag(obj.mean_raw(1)));
 %% ���� 2
 temp = obj.raw(obj.zone2(end)-m:obj.zone2(end),:);
 %������� ������ ������� �������� (��������� ��� ���������)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %���������� ���. �������� �� ������� (1)
 obj.mean_raw(2) = mean(temp(:));
 %���������� ��������� �� ������� (2)-(3)
 obj.disp_re(2) = var(real(temp(:))); %��������� ������������ ����� �����
 obj.disp_im(2) = var(imag(temp(:))); %��������� ������ ����� �����
 %����������� ��������� �� ������� (4)
 obj.raw(obj.zone2,:) = 1/sqrt(obj.disp_re(2))*(real(obj.raw(obj.zone2,:))-real(obj.mean_raw(2)))
+...
 1j/sqrt(obj.disp_im(2))*(imag(obj.raw(obj.zone2,:))-imag(obj.mean_raw(2)));
 %% ���� 3
 temp = obj.raw(obj.zone3(end)-m:obj.zone3(end),:);
 %������� ������ ������� �������� (��������� ��� ���������)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %���������� ���. �������� �� ������� (1)
 obj.mean_raw(3) = mean(temp(:));
 %���������� ��������� �� ������� (2)-(3)
 obj.disp_re(3) = var(real(temp(:))); %��������� ������������ ����� �����
 obj.disp_im(3) = var(imag(temp(:))); %��������� ������ ����� �����
 %����������� ��������� �� ������� (4)
 obj.raw(obj.zone3,:) = 1/sqrt(obj.disp_re(3))*(real(obj.raw(obj.zone3,:))-real(obj.mean_raw(3)))
+...
 1j/sqrt(obj.disp_im(3))*(imag(obj.raw(obj.zone3,:))-imag(obj.mean_raw(3)));
 %������������� ���������
 obj.flag_raf = 1;
 end
 end
end
