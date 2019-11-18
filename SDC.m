classdef SDC < handle
    %SDC осуществляет селекцию движущищхся целей

    properties
 
        m_ConvData %указатель на объект класса ConvData

 
        %определяем параметры когерентной СДЦ
 
        n_chain = 8 %число отсчетов в пачке
 
        n_step = 4 %шаг смещения пачки
 
        wnd %весовое окно на пачку
 
        zCh = -2:2; %номера нулевых каналов СДЦ (содержат отражания от земли и стац. объектов)

 
        sdc_data %выходной массив данных
 
        sdc_ang %выходной массив кодов углов данных

 
        %параметры отображения сигнала
 
        med %медиана амплитуд данных на входе
 
    end

    methods
        function obj = SDC(hConvData) %конструктор объекта
       
            obj.m_ConvData = hConvData; %создаем связь с объектом класса ConvData

 
            %создаем окно для наложения на пачку
 
            obj.wnd = ones(1,obj.n_chain);

 
            %определяем параметры отображения
 
            obj.med = median(abs(obj.m_ConvData.convolved_data(:)));
 
        end
        
        function chpk1(obj) %однократное череспериодное вычитание
 
            obj.sdc_data = []; %удаляем результаты прошлых вычислений
 
            temp = obj.m_ConvData.convolved_data; %исходные данны для сдц
 
            n = size(temp,2); %число периодов излучения сигнала
 
            obj.sdc_data = abs(temp(:,2:n)) - abs(temp(:,1:n-1));
 
            obj.sdc_ang = 2:n;

 
            obj.show;
 
        end
        
        function sdc_fft(obj) %когерентнаф сдц через цифровую фильтрацию
 
            obj.sdc_data = []; %удаляем результаты прошлых вычислений
 
            temp = obj.m_ConvData.convolved_data; %исходные данны для свертки
 
            n = size(temp,2); %число периодов по дальности
 
            n = fix(n/obj.n_step)*obj.n_step; %общее число обрабатываемых каналов
 
            temp(:,(n+1):end) = []; %удаляем на обрабатываемые данные
 
            obj.sdc_data = zeros(size(temp,1),obj.n_chain,n/obj.n_step-1);
 
            for ii = 1:n/obj.n_step-1
 
                temp2 = fft(temp(:,(ii-1)*obj.n_step+(1:obj.n_chain)),[],2); %формула (12)
 
                obj.sdc_data(:,:,ii) = fftshift(temp2,2);
 
            end
            
            obj.sdc_ang = 1:obj.n_step:n-obj.n_step;

 
            obj.show;
 
        end
        
        function sdc_mo(obj)
            
            obj.sdc_data = []; % удаляем результаты прошлых вычислений
            
            temp = obj.m_ConvData.convolved_data; % исходные данные для свертки
            
            n = size(temp, 2); % xb
            
            obj.sdc_data = abs(temp(:,2:n)) - abs(obj.form_map_data);
            
            obj.show;
       
        end
        
        function show(obj)
 
            figure;
 
            subplot(1,2,1);
 
            imagesc(abs(obj.m_ConvData.convolved_data),[0 10*abs(obj.med)]);
 
            xlabel('Угол, отс');
 
            ylabel('Дальность, отс');
 
            subplot(1,2,2);
 
            %проверка типа СДЦ
 
            sz = size(obj.sdc_data);
 
            if length(sz)>2
 
                %когерентная сдц
 
                viewdata = obj.sdc_data;
 
                viewdata(:,obj.n_chain/2+obj.zCh,:) = 0; %формула (13)
 
                viewdata = ifft(ifftshift(viewdata,2),[],2); %формула (14)
 
                viewdata = viewdata(:,obj.n_chain/2,:);
 
                viewdata = permute(viewdata,[1,3,2]);
 
                imagesc(abs(viewdata),[0 10*abs(obj.med)]);
 
            else
                
                %некогерентая сдц
 
                imagesc(abs(obj.sdc_data),[0 10*abs(obj.med)]);
 
            end
            
            xlabel('Угол, отс');
 
            ylabel('Дальность, отс');
 
        end
        
        function Kp = calc_kp(obj,r,a)
 
            %оценка коэффициента подавления в выбранной точке
 
            %проверка типа СДЦ
 
            sz = size(obj.sdc_data);
 
            if length(sz)>2
 
                viewdata = obj.sdc_data;
 
                viewdata(:,obj.n_chain/2+obj.zCh,:) = 0; %формула (13)
 
                viewdata = ifft(ifftshift(viewdata,2),[],2); %формула (14)
 
                viewdata = viewdata(:,obj.n_chain/2,:);
 
                viewdata = permute(viewdata, [1, 3, 2]);
 
                amp_sdc = interp1(obj.sdc_ang, viewdata(r,:), a);
 
                Kp = 20*log10(abs(amp_sdc)/abs(obj.m_ConvData.convolved_data(r,a)));
 
            else
                
                Kp = 20*log10(abs(obj.sdc_data(r,a))/abs(obj.m_ConvData.convolved_data(r,a)));
 
            end
            
        end
        
        function map = form_map_data(obj)
 
            %формирование исходных данных для карты местности по сдц
 
            map = obj.sdc_data;
 
            idx = ones(1,obj.n_chain);
 
            idx(obj.n_chain/2+obj.zCh) = 0;
 
            map(:,boolean(idx),:) = 0; %формула (13)
 
            map = ifft(ifftshift(map,2),[],2); %формула (14)
 
            map = map(:,obj.n_chain/2,:);
 
            map = permute(map,[1,3,2]);
 
        end
        
    end
end