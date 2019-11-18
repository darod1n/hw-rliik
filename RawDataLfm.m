classdef RawDataLfm < handle           %наследник глобального класса hanlde в матлабе
                                    %RawData описывает загрузка данных РЛС и их "рафинирование"
    properties
        fname %имя файла
        fpath %путь к файлу
        raw %массив амплитуд в координатах дальность-угол
        angles %массив кодов углов в отсчетах
        m %медиана абсолютных амплитуд в файле
        flag_raf%флаг рафинирования
        mean_raw = zeros(1,3); %мат ожидание шумов
        disp_re = zeros(1,3); %дисперсия вещественной части шумов
        disp_im = zeros(1,3); %дисперсия мнимой части шумов

        %выделяем зоны для обработки
        zone1 = [1:1000];
        zone2 = [1001:2000];
        zone3 = [2001:16000];
    end

    methods
        function obj = RawDataLfm(fpath,fn) %конструктор класса
            if nargin < 1
                obj.open;
            else
                obj.open(fpath,fn);
            end
            %отображение данных для того, чтоб убедиться, что в области
            %больших дальностей присутствуют только шумы
            obj.show;
        end
        function open(obj,fpath,fn) %открытие файла для обработки
            if nargin < 1
            %выбор файла для обработки
 [obj.fname, obj.fpath] = uigetfile('*.mat');
 else
 obj.fpath = fpath;
 obj.fname = fn;
 end
 %приведение типов к float
 temp = (load([obj.fpath obj.fname],'lfmraw'));
 obj.raw = single(temp.lfmraw);
 temp = (load([obj.fpath obj.fname],'lfmangles'));
 obj.angles = single(temp.lfmangles);
 %вычисление медианы для отображения данных
 obj.m = median(abs(obj.raw(:)));
 %сброс флага рафинирования
 obj.flag_raf = 0; %обработка не проводилась
 end
 function show(obj) %отображение данных
 figure;
 imagesc(abs(obj.raw),[0 2*obj.m]);
 end
 function rafinate(obj,m)
 %выборка данных: все каналы дальности и m последних каналов
 %дальности
 %% зона 1
 temp = obj.raw(obj.zone1(end)-m:obj.zone1(end),:);
 %выборка только шумовых отсчетов (исключаем все нешумовые)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %вычисление мат. ожидания по формуле (1)
 obj.mean_raw(1) = mean(temp);
 %вычисление дисперсий по формуле (2)-(3)
 obj.disp_re(1) = var(real(temp)); %дисперсия вещественной части шумов
 obj.disp_im(1) = var(imag(temp)); %дисперсия мнимой части шумов
 %компенсация квадратур по фомруле (4)
 obj.raw(obj.zone1,:) = 1/sqrt(obj.disp_re(1))*(real(obj.raw(obj.zone1,:))-real(obj.mean_raw(1)))
+...
 1j/sqrt(obj.disp_im(1))*(imag(obj.raw(obj.zone1,:))-imag(obj.mean_raw(1)));
 %% зона 2
 temp = obj.raw(obj.zone2(end)-m:obj.zone2(end),:);
 %выборка только шумовых отсчетов (исключаем все нешумовые)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %вычисление мат. ожидания по формуле (1)
 obj.mean_raw(2) = mean(temp(:));
 %вычисление дисперсий по формуле (2)-(3)
 obj.disp_re(2) = var(real(temp(:))); %дисперсия вещественной части шумов
 obj.disp_im(2) = var(imag(temp(:))); %дисперсия мнимой части шумов
 %компенсация квадратур по фомруле (4)
 obj.raw(obj.zone2,:) = 1/sqrt(obj.disp_re(2))*(real(obj.raw(obj.zone2,:))-real(obj.mean_raw(2)))
+...
 1j/sqrt(obj.disp_im(2))*(imag(obj.raw(obj.zone2,:))-imag(obj.mean_raw(2)));
 %% зона 3
 temp = obj.raw(obj.zone3(end)-m:obj.zone3(end),:);
 %выборка только шумовых отсчетов (исключаем все нешумовые)
 temp = temp(temp(:) < median(temp(:)));
 temp = temp(:);
 %вычисление мат. ожидания по формуле (1)
 obj.mean_raw(3) = mean(temp(:));
 %вычисление дисперсий по формуле (2)-(3)
 obj.disp_re(3) = var(real(temp(:))); %дисперсия вещественной части шумов
 obj.disp_im(3) = var(imag(temp(:))); %дисперсия мнимой части шумов
 %компенсация квадратур по фомруле (4)
 obj.raw(obj.zone3,:) = 1/sqrt(obj.disp_re(3))*(real(obj.raw(obj.zone3,:))-real(obj.mean_raw(3)))
+...
 1j/sqrt(obj.disp_im(3))*(imag(obj.raw(obj.zone3,:))-imag(obj.mean_raw(3)));
 %рафинирование выполнено
 obj.flag_raf = 1;
 end
 end
end
