classdef SDC < handle
    %SDC ������������ �������� ����������� �����

    properties
 
        m_ConvData %��������� �� ������ ������ ConvData

 
        %���������� ��������� ����������� ���
 
        n_chain = 8 %����� �������� � �����
 
        n_step = 4 %��� �������� �����
 
        wnd %������� ���� �� �����
 
        zCh = -2:2; %������ ������� ������� ��� (�������� ��������� �� ����� � ����. ��������)

 
        sdc_data %�������� ������ ������
 
        sdc_ang %�������� ������ ����� ����� ������

 
        %��������� ����������� �������
 
        med %������� �������� ������ �� �����
 
    end

    methods
        function obj = SDC(hConvData) %����������� �������
       
            obj.m_ConvData = hConvData; %������� ����� � �������� ������ ConvData

 
            %������� ���� ��� ��������� �� �����
 
            obj.wnd = ones(1,obj.n_chain);

 
            %���������� ��������� �����������
 
            obj.med = median(abs(obj.m_ConvData.convolved_data(:)));
 
        end
        
        function chpk1(obj) %����������� �������������� ���������
 
            obj.sdc_data = []; %������� ���������� ������� ����������
 
            temp = obj.m_ConvData.convolved_data; %�������� ����� ��� ���
 
            n = size(temp,2); %����� �������� ��������� �������
 
            obj.sdc_data = abs(temp(:,2:n)) - abs(temp(:,1:n-1));
 
            obj.sdc_ang = 2:n;

 
            obj.show;
 
        end
        
        function sdc_fft(obj) %����������� ��� ����� �������� ����������
 
            obj.sdc_data = []; %������� ���������� ������� ����������
 
            temp = obj.m_ConvData.convolved_data; %�������� ����� ��� �������
 
            n = size(temp,2); %����� �������� �� ���������
 
            n = fix(n/obj.n_step)*obj.n_step; %����� ����� �������������� �������
 
            temp(:,(n+1):end) = []; %������� �� �������������� ������
 
            obj.sdc_data = zeros(size(temp,1),obj.n_chain,n/obj.n_step-1);
 
            for ii = 1:n/obj.n_step-1
 
                temp2 = fft(temp(:,(ii-1)*obj.n_step+(1:obj.n_chain)),[],2); %������� (12)
 
                obj.sdc_data(:,:,ii) = fftshift(temp2,2);
 
            end
            
            obj.sdc_ang = 1:obj.n_step:n-obj.n_step;

 
            obj.show;
 
        end
        
        function sdc_mo(obj)
            
            obj.sdc_data = []; % ������� ���������� ������� ����������
            
            temp = obj.m_ConvData.convolved_data; % �������� ������ ��� �������
            
            n = size(temp, 2); % xb
            
            obj.sdc_data = abs(temp(:,2:n)) - abs(obj.form_map_data);
            
            obj.show;
       
        end
        
        function show(obj)
 
            figure;
 
            subplot(1,2,1);
 
            imagesc(abs(obj.m_ConvData.convolved_data),[0 10*abs(obj.med)]);
 
            xlabel('����, ���');
 
            ylabel('���������, ���');
 
            subplot(1,2,2);
 
            %�������� ���� ���
 
            sz = size(obj.sdc_data);
 
            if length(sz)>2
 
                %����������� ���
 
                viewdata = obj.sdc_data;
 
                viewdata(:,obj.n_chain/2+obj.zCh,:) = 0; %������� (13)
 
                viewdata = ifft(ifftshift(viewdata,2),[],2); %������� (14)
 
                viewdata = viewdata(:,obj.n_chain/2,:);
 
                viewdata = permute(viewdata,[1,3,2]);
 
                imagesc(abs(viewdata),[0 10*abs(obj.med)]);
 
            else
                
                %������������ ���
 
                imagesc(abs(obj.sdc_data),[0 10*abs(obj.med)]);
 
            end
            
            xlabel('����, ���');
 
            ylabel('���������, ���');
 
        end
        
        function Kp = calc_kp(obj,r,a)
 
            %������ ������������ ���������� � ��������� �����
 
            %�������� ���� ���
 
            sz = size(obj.sdc_data);
 
            if length(sz)>2
 
                viewdata = obj.sdc_data;
 
                viewdata(:,obj.n_chain/2+obj.zCh,:) = 0; %������� (13)
 
                viewdata = ifft(ifftshift(viewdata,2),[],2); %������� (14)
 
                viewdata = viewdata(:,obj.n_chain/2,:);
 
                viewdata = permute(viewdata, [1, 3, 2]);
 
                amp_sdc = interp1(obj.sdc_ang, viewdata(r,:), a);
 
                Kp = 20*log10(abs(amp_sdc)/abs(obj.m_ConvData.convolved_data(r,a)));
 
            else
                
                Kp = 20*log10(abs(obj.sdc_data(r,a))/abs(obj.m_ConvData.convolved_data(r,a)));
 
            end
            
        end
        
        function map = form_map_data(obj)
 
            %������������ �������� ������ ��� ����� ��������� �� ���
 
            map = obj.sdc_data;
 
            idx = ones(1,obj.n_chain);
 
            idx(obj.n_chain/2+obj.zCh) = 0;
 
            map(:,boolean(idx),:) = 0; %������� (13)
 
            map = ifft(ifftshift(map,2),[],2); %������� (14)
 
            map = map(:,obj.n_chain/2,:);
 
            map = permute(map,[1,3,2]);
 
        end
        
    end
end