%% @Time:2025/11/02
%% @Matlab code for the paper《Simultaneous Calibration of Multicoordinates for a Dual-Robot System by Solving  the AXB = YCZ Problem》.
%% @Rewritten from the Julia code of the GC_DualRobot project.
%% Reproduced author: Wuyi, Phd student of Hebut. email:wuyiii2022@163.com 
%% mian script

% 读取数据
fprintf('步骤1: 读取输入文件...\n');
[A, B, C, n_valid] = read_and_convert_valid_data('A_1.txt','B_1.txt','C_1.txt');

fprintf('成功读取 %d 组有效数据\n', n_valid);
if n_valid < 3
    error('有效数据不足，至少需要3组数据才能进行校准');
end

% 迭代参数设置
% confₘ₂ = ((max_iter = 200, m = 3, η=0.1,μ=0.1,scale=1, err = "NAN", stop = true, reg = true, τ=0.3))	

conf.max_iter = 200;
conf.m =3;
conf.eta = 0.1;
conf.scale = 1;
conf.stop = true;
conf.err = "LS";
conf.reg = true;
conf.tau = 0.3;
conf.mu = [1,1e-5,1,1,1];
conf.svd = false;

% 闭式解法
[X_init, Y_init, Z_init] = AXB_YCZ_Close(A, B, C, conf);

% 迭代解法
V = {X_init,Y_init, Z_init}; P = {A,B,C};
[Result,err] = iter_solve_Katyusha(V,P,conf);

% 检查配置中是否启用SVD选项
    if isfield(conf, 'svd') && conf.svd
        % 获取数据点数量
        n = length(A);
        
        % 从变换矩阵中提取旋转部分（3x3）和平移部分（3x1）
        RA = cell(1, n);
        RB = cell(1, n);
        RC = cell(1, n);
        tA = cell(1, n);
        tB = cell(1, n);
        tC = cell(1, n);
        
        for i = 1:n
            % 提取旋转矩阵部分（左上角3x3）
            RA{i} = A{i}(1:3, 1:3);
            RB{i} = B{i}(1:3, 1:3);
            RC{i} = C{i}(1:3, 1:3);
            
            % 提取平移向量部分（第4列的前3个元素）
            tA{i} = A{i}(1:3, 4);
            tB{i} = B{i}(1:3, 4);
            tC{i} = C{i}(1:3, 4);
        end
        
        % 从优化结果中提取旋转矩阵并正交化
        % Result是3x12矩阵，结构为 [Rx, Ry, Rz, Tx, Ty, Tz]
        Rx_sln = Result(1:3, 1:3);      % 提取Rx部分
        Ry_sln = Result(1:3, 4:6);      % 提取Ry部分
        Rz_sln = Result(1:3, 7:9);      % 提取Rz部分
        
        % 投影到SO(3)流形（确保旋转矩阵正交且行列式为+1）
        Rx_sln_ortho = ProjectToSO3(Rx_sln);
        Ry_sln_ortho = ProjectToSO3(Ry_sln);
        Rz_sln_ortho = ProjectToSO3(Rz_sln);
        
        % 使用SVD方法求解平移向量
        % 假设T_SVD函数已实现，返回三个变换矩阵
        [X_r, Y_r, Z_r] = T_SVD(Rx_sln_ortho, Ry_sln_ortho, Rz_sln_ortho, RA, RB, RC, tA, tB, tC);
        
        % 显示结果
        display(X_r);
        display(Y_r);
        display(Z_r);
%         display(err)
    end
    
    % 如果未启用SVD，使用默认后处理方法
    [X_r, Y_r, Z_r] = post_process_Wang(Result);
    % 显示结果
    display(X_r);
    display(Y_r);
    display(Z_r);
%     display(err)
    figure;
    plot(err,'LineWidth',1.5);
    xlabel('Number of iterations');
    ylabel('Mean error of ||AXB-YCZ||_2');

%% ----------------------- 子函数 ----------------------------------------- %%
function T = eulerZYXToRotationMatrix(eulerAngles)
% eulerZYXToRotationMatrix - 将ZYX欧拉角转换为旋转矩阵
% 输入:
%   eulerAngles - 1x6或6x1向量，包含[X,Y,Z,R,P,Y]（角度）                
% 输出:
%   T - 4x4旋转矩阵

    % 输入验证
    if numel(eulerAngles) ~= 6
        error('输入必须包含6个元素');
    end
    
    % 确保输入为行向量
    eulerAngles = eulerAngles(:)';  % 转换为行向量
    
    % 提取各角度（ZYX顺序）
    gamma = eulerAngles(6)/180*pi;  % Z轴旋转角
    beta  = eulerAngles(5)/180*pi;  % Y轴旋转角
    alpha = eulerAngles(4)/180*pi;  % X轴旋转角
    
    % 计算三角函数值（预计算提高效率）
    cos_g = cos(gamma); sin_g = sin(gamma);
    cos_b = cos(beta);  sin_b = sin(beta);
    cos_a = cos(alpha); sin_a = sin(alpha);
    
    % 计算绕各轴的基本旋转矩阵
    % 绕Z轴旋转矩阵
    Rz = [cos_g, -sin_g, 0;
          sin_g,  cos_g, 0;
          0,      0,     1];
    
    % 绕Y轴旋转矩阵
    Ry = [cos_b,  0, sin_b;
          0,      1, 0;
          -sin_b, 0, cos_b];
    
    % 绕X轴旋转矩阵
    Rx = [1, 0,      0;
          0, cos_a, -sin_a;
          0, sin_a,  cos_a];
    
    % 组合旋转矩阵（ZYX顺序：先Z后Y最后X）
    R = Rz * Ry * Rx;
    
    % 数值稳定性检查（可选）
    if abs(det(R) - 1) > 1e-10
        warning('旋转矩阵行列式偏离1较多: %.6e', det(R));
    end
    
    t = [eulerAngles(1),eulerAngles(2),eulerAngles(3)]';
    % 合并齐次矩阵
    T = [R,t;0,0,0,1];
end

%% 数据读取函数（保持不变）
function [A_cell, B_cell, C_cell, n_valid] = read_and_convert_valid_data(a,b,c)
    % 读取有效数据文件
    A_data = dlmread(a);
    B_data = dlmread(b);
    C_data = dlmread(c);
    
    n_valid = size(A_data, 1);
    
    % 转换为齐次变换矩阵
    A_cell = cell(1, n_valid);
    B_cell = cell(1, n_valid);
    C_cell = cell(1, n_valid);
    
    for i = 1:n_valid
        A_cell{i} = eulerZYXToRotationMatrix(A_data(i, :));
        A_cell{i}(1:3,4)=A_cell{i}(1:3,4)./1000;
        B_cell{i} = reshape(B_data(i, :), 4, 4)';
        B_cell{i}(1:3,4)=B_cell{i}(1:3,4);
        C_cell{i} = eulerZYXToRotationMatrix(C_data(i, :));
        C_cell{i}(1:3,4)=C_cell{i}(1:3,4)./1000;
    end
end

%% 
function res = toSM(M)
  %%Julia 中将矩阵转换为静态矩阵
  res = [M(1,1) M(1,2) M(1,3) M(1,4); M(2,1) M(2,2) M(2,3) M(2,4); M(3,1) M(3,2) M(3,3) M(3,4); 0.0 0.0 0.0 1.0];
end
%%
function res = toV(Init)
   % 将三个齐次变换矩阵转换为优化变量矩阵V
   % 输入: Init - 元胞数组，包含三个4x4变换矩阵 {X, Y, Z}
   % 输出: V - 3x12矩阵，包含旋转和平移部分
    
   % 提取旋转矩阵部分（3x3）
   R1 = Init{1}(1:3, 1:3);
   R2 = Init{2}(1:3, 1:3);
   R3 = Init{3}(1:3, 1:3);
   
   % 提取平移向量部分（3x1）
   t1 = Init{1}(1:3, 4);
   t2 = Init{2}(1:3, 4);
   t3 = Init{3}(1:3, 4);
   
   % 水平连接所有部分
   res = [R1, R2, R3, t1, t2, t3]; 
end
%%
function R = projectToSO3(M)
    % 将输入矩阵正交化
    % 输入: M ，3x3旋转矩阵
    % 输出: R ，正交化之后的3x3旋转矩阵
    [U, ~, V] = svd(M);
    R = U * V';
    if det(R) < 0
        V(:,3) = -V(:,3);
        R = U * V';
    end
end
%%
function V_out = normfunc(V)
    % 将优化变量矩阵V中的三个旋转矩阵块进行正交化处理
    % 输入: V ，3x12 拼接矩阵[R_X, R_Y, R_Z, t_X, t_Y, t_Z]
    % 输出: V_out ，正交化之后的3x3旋转矩阵，
    V_out = V;
	V_out(1:3,1:3) = projectToSO3(V(1:3,1:3));
	V_out(1:3,4:6) = projectToSO3(V(1:3,4:6));
	V_out(1:3,7:9) = projectToSO3(V(1:3,7:9));
end
%%
function v = tovec(R)
    v = R(:);  % MATLAB中的冒号操作符将矩阵转换为列向量
end
%%
function Mc = toMc(vecRc)
    Mc = zeros(9, 81);
    for i = 1:9
        startCol = (i-1)*9 + 1;
        endCol = i*9;
        Mc(i, startCol:endCol) = vecRc(:).'; % 确保行向量赋值
    end
end

%%
function R = solveR(A, B, C)
% 闭式解法求解AXB=YCZ问题中的旋转矩阵
% 输入:
%   A, B, C - 元胞数组，每个元素为3x3旋转矩阵
%             A{i}: 机器人1的基座到末端法兰变换的旋转部分
%             B{i}: 工具或传感器的旋转矩阵
%             C{i}: 机器人2的末端法兰到工具的旋转矩阵
% 输出:
%   R - 3x3旋转矩阵

    % 获取数据点数量
    n = length(A);
    
    % 验证输入一致性
    if length(A) ~= length(B) || length(A) ~= length(C)
        error('输入数据A、B、C的长度必须一致');
    end
    
    % 初始化大矩阵M（9n x 90）
    M = zeros(9*n, 90);
    
    % 构造线性系统矩阵M
    for i = 1:n
        % 计算当前数据点的块矩阵
        kron_part = kron(B{i}', A{i});  % Kronecker积: B^T ⊗ A
        mc_part = toMc(tovec(C{i}));    % 矩阵C的块对角形式
        
        % 组合当前数据点的贡献矩阵（9x90）
        block_matrix = [kron_part, -mc_part];
        
        % 将块矩阵放入M的相应行
        start_row = 9*(i-1) + 1;
        end_row = 9*i;
        M(start_row:end_row, :) = block_matrix;
    end
    
    % SVD分解求解最小二乘问题
    [~, ~, V] = svd(M' * M);
    
    % 提取最小奇异值对应的特征向量并重塑为旋转矩阵
    v_min = V(1:9, 90);                    % 取第90列的前9个元素
    v_min_normalized = v_min / norm(v_min); % 归一化
    R_reshaped = reshape(v_min_normalized, 3, 3); % 重塑为3x3矩阵
    R = R_reshaped * sqrt(3);               % 应用缩放因子√3
end

%%
function [R1, R2, R3] = solveRS(RA, RB, RC)
% 闭式解法求解三个旋转矩阵
% 输入: RA, RB, RC - 元胞数组，每个元素为3x3旋转矩阵
% 输出: R1, R2, R3 - 3x3旋转矩阵估计

    % 调用solveR函数三次，使用不同输入组合
    R1 = solveR(RA, RB, RC);  % 原始输入
    
    % 对RA的每个元素转置（类似Julia的transpose.操作）
    RA_transposed = cellfun(@transpose, RA, 'UniformOutput', false);
    R2 = solveR(RA_transposed, RC, RB);  % 输入顺序：转置RA、RC、RB
    
    % 对RB的每个元素转置
    RB_transposed = cellfun(@transpose, RB, 'UniformOutput', false);
    R3 = solveR(RC, RB_transposed, RA);  % 输入顺序：RC、转置RB、RA

end

%%
function [Tx, Ty, Tz, err] = solveTs(RA, TA, RB, TB, RC, TC, Rx, Ry, Rz)
% 求解平移向量和计算校准误差
% 输入:
%   RA, RB, RC - 元胞数组，每个元素为3x3旋转矩阵
%   TA, TB, TC - 元胞数组，每个元素为3x1平移向量
%   Rx, Ry, Rz - 3x3矩阵，旋转矩阵估计值
% 输出:
%   Tx, Ty, Tz - 3x1平移向量
%   err - 标量，校准误差

    % 获取数据点数量
    n = length(RA);
    
    % 验证输入一致性
    if length(RA) ~= length(TA) || length(RA) ~= length(RB) || ...
       length(RA) ~= length(TB) || length(RA) ~= length(RC) || length(RA) ~= length(TC)
        error('输入数组长度必须一致');
    end
    
    % 初始化雅可比矩阵J和向量b
    J = zeros(3*n, 9);
    b = zeros(3*n, 1);
    
    % 构造线性系统
    for i = 1:n
        % 提取当前数据点
        RA_i = RA{i};
        TA_i = TA{i};
        RB_i = RB{i};
        TB_i = TB{i};
        RC_i = RC{i};
        TC_i = TC{i};
        
        % 填充雅可比矩阵J（3x9块）
        start_row = 3*i - 2;
        end_row = 3*i;
        J_block = [RA_i, -eye(3), -Ry * RC_i];
        J(start_row:end_row, :) = J_block;
        
        % 填充向量b（3x1块）
        b_block = Ry * TC_i - TA_i - RA_i * Rx * TB_i;
        b(start_row:end_row, 1) = b_block;
    end
    
    % 最小二乘求解（使用正规方程）
    t = (J' * J) \ (J' * b);
    
    % 提取平移向量
    Tx = t(1:3);
    Ty = t(4:6);
    Tz = t(7:9);
    
    % 使用第一个数据点构造齐次变换矩阵（用于误差计算）
    A_p = [RA{1}, TA{1}; 0, 0, 0, 1];
    B_p = [RB{1}, TB{1}; 0, 0, 0, 1];
    C_p = [RC{1}, TC{1}; 0, 0, 0, 1];
    X_p = [Rx, Tx; 0, 0, 0, 1];
    Y_p = [Ry, Ty; 0, 0, 0, 1];
    Z_p = [Rz, Tz; 0, 0, 0, 1];
    
    % 计算误差（2-范数）
    err = norm(A_p * X_p * B_p - Y_p * C_p * Z_p, 2);
end

%%
function [Rx, Tx, Ry, Ty, Rz, Tz] = solveBestT(RA, TA, RB, TB, RC, TC, Rx_initial, Ry_initial, Rz_initial)
    % 调整旋转矩阵的符号以确保行列式为+1
    Rx = sign(det(Rx_initial)) * Rx_initial; %Rx_initial;%
    Ry = sign(det(Ry_initial)) * Ry_initial; %Ry_initial;%
    Rz = sign(det(Rz_initial)) * Rz_initial; %Rz_initial;%

    % 调用solveTs函数求解平移向量（假设已有solveTs实现）
    [Tx, Ty, Tz, err] = solveTs(RA, TA, RB, TB, RC, TC, Rx, Ry, Rz);
    display(err);
end

%%
function T = rt2T(R, t)
% 将旋转矩阵和平移向量组合为齐次变换矩阵
% 输入: R - 3x3旋转矩阵, t - 3x1平移向量
% 输出: T - 4x4齐次变换矩阵
    T = [R, t; 0, 0, 0, 1];
end

%%
function [X, Y, Z] = AXB_YCZ_Close(A, B, C, conf)
% AXB=YCZ校准问题的闭式解法
% 输入:
%   A, B, C - 元胞数组，每个元素为4x4齐次变换矩阵
%   conf - 配置结构体（可选，用于扩展参数）
% 输出:
%   X, Y, Z - 4x4齐次变换矩阵，校准结果

    % 获取数据点数量
    n = length(A);
    
    % 验证输入一致性
    if length(A) ~= length(B) || length(A) ~= length(C)
        error('输入数据A、B、C的长度必须一致');
    end
    
    % 从齐次变换矩阵中提取旋转矩阵部分（3x3）
    RA = cell(1, n);
    RB = cell(1, n);
    RC = cell(1, n);
    for i = 1:n
        RA{i} = A{i}(1:3, 1:3);
        RB{i} = B{i}(1:3, 1:3);
        RC{i} = C{i}(1:3, 1:3);
    end
    
    % 从齐次变换矩阵中提取平移向量部分（3x1）
    tA = cell(1, n);
    tB = cell(1, n);
    tC = cell(1, n);
    for i = 1:n
        tA{i} = A{i}(1:3, 4);
        tB{i} = B{i}(1:3, 4);
        tC{i} = C{i}(1:3, 4);
    end
    
    % 调用solveRS函数求解旋转矩阵初始估计
    [Rx_init, Ry_init, Rz_init] = solveRS(RA, RB, RC);
    
%     [X, Y, Z] = T_SVD(Rx_init, Ry_init, Rz_init, RA, RB, RC, tA, tB, tC);  % SVD方法求解平移矢量


    % 调用solveBestT函数求解平移向量并处理符号歧义
     [Rx, Tx, Ry, Ty, Rz, Tz] = solveBestT(RA, tA, RB, tB, RC, tC, Rx_init, Ry_init, Rz_init);
    
    % 姿态正交处理
%     Rx = projectToSO3(Rx);
%     Ry = projectToSO3(Ry);
%     Rz = projectToSO3(Rz);
    % 将旋转和平移组合为齐次变换矩阵
     X = rt2T(Rx, Tx);
     Y = rt2T(Ry, Ty);
     Z = rt2T(Rz, Tz);
end

%%
function J = GetJ(RA_noise, RY_sln, RC_noise)
% GetJ - 构造雅可比矩阵J用于平移向量求解
% 输入:
%   RA_noise, RC_noise - 元胞数组，每个元素为3x3旋转矩阵（含噪声的测量数据）
%   RY_sln - 3x3旋转矩阵（已校准的Y旋转矩阵）
% 输出:
%   J - 3M×9雅可比矩阵，M为数据点数量

    M = length(RA_noise);
    J = zeros(3*M, 9);
    
    for i = 1:M
        % 构造当前数据点的雅可比块
        J_block = [RA_noise{i}, -eye(3), -RY_sln * RC_noise{i}];
        
        % 将块放入雅可比矩阵的相应位置
        start_row = 3*i - 2;
        end_row = 3*i;
        J(start_row:end_row, :) = J_block;
    end
end

function p = Getp(RA_noise, RX_sln, RY_sln, tA_noise, tB_noise, tC_noise)
% Getp - 构造残差向量p用于平移向量求解
% 输入:
%   RA_noise - 元胞数组，每个元素为3x3旋转矩阵
%   RX_sln, RY_sln - 3x3旋转矩阵（已校准的旋转矩阵）
%   tA_noise, tB_noise, tC_noise - 元胞数组，每个元素为3x1平移向量
% 输出:
%   p - 3M×1残差向量，M为数据点数量

    M = length(RA_noise);
    p = zeros(3*M, 1);
    
    for i = 1:M
        % 计算当前数据点的残差
        p_i = -tA_noise{i} - RA_noise{i} * RX_sln * tB_noise{i} + RY_sln * tC_noise{i};
        
        % 将残差放入向量的相应位置
        start_row = 3*i - 2;
        end_row = 3*i;
        p(start_row:end_row) = p_i;
    end
end

function [X_s, Y_s, Z_s] = T_SVD(RX_sln, RY_sln, RZ_sln, RA, RB, RC, tA, tB, tC)
% T_SVD - 通过SVD方法求解平移向量并构造齐次变换矩阵
% 输入:
%   RX_sln, RY_sln, RZ_sln - 3x3旋转矩阵（已求解的旋转部分）
%   RA, RB, RC - 元胞数组，每个元素为3x3旋转矩阵（测量数据）
%   tA, tB, tC - 元胞数组，每个元素为3x1平移向量（测量数据）
% 输出:
%   X_s, Y_s, Z_s - 4x4齐次变换矩阵（完整的校准结果）

    % 获取数据点数量
    M = length(RA);
    
    % 验证输入一致性
    if length(RA) ~= length(RB) || length(RA) ~= length(RC) || ...
       length(RA) ~= length(tA) || length(RA) ~= length(tB) || length(RA) ~= length(tC)
        error('所有输入数组的长度必须一致');
    end
    
    % 构造雅可比矩阵J
    J = GetJ(RA, RY_sln, RC);
    
    % 构造残差向量p
    p = Getp(RA, RX_sln, RY_sln, tA, tB, tC);
    
    % 最小二乘求解平移向量
    t_sln = J \ p;
    
    % 提取各个平移分量
    tX_sln = t_sln(1:3);
    tY_sln = t_sln(4:6);
    tZ_sln = t_sln(7:9);
    
    % 构造齐次变换矩阵
    X_s = [RX_sln, tX_sln; 0, 0, 0, 1];
    Y_s = [RY_sln, tY_sln; 0, 0, 0, 1];
    Z_s = [RZ_sln, tZ_sln; 0, 0, 0, 1];
end

%%
function error_val = errN(V, P, grad_V, eta, conf)
% errN - 空误差函数（用于测试或忽略误差计算）
% 输入: V - 优化变量, P - 数据, grad_V - 梯度, eta - 步长, conf - 配置
% 输出: error_val - 固定返回0
    error_val = 0.0;
end

%%
function [X, Y, Z] = post_process_Wang(V_t)
% post_process_Wang - 将优化变量向量解码为齐次变换矩阵
% 输入:
%   V_t - 3x12矩阵，优化变量向量 [Rx,Ry,Rz,Tx,Ty,Tz]
%         第1-3列: Rx (3x3旋转矩阵的向量化)
%         第4-6列: Ry (3x3旋转矩阵的向量化)  
%         第7-9列: Rz (3x3旋转矩阵的向量化)
%         第10-12列: Tx,Ty,Tz (3x1平移向量)
% 输出:
%   X, Y, Z - 4x4齐次变换矩阵

    % 输入验证
    if size(V_t, 1) ~= 3 || size(V_t, 2) ~= 12
        error('输入矩阵V_t必须是3x12维度');
    end
    
    % 提取旋转矩阵部分
    Rx = V_t(:, 1:3);   % 前3列对应Rx
    Ry = V_t(:, 4:6);   % 中间3列对应Ry
    Rz = V_t(:, 7:9);   % 后3列对应Rz
    
    % 提取平移向量部分
    Tx = V_t(:, 10);
    Ty = V_t(:, 11);
    Tz = V_t(:, 12);
    
    % 对旋转矩阵进行正交化处理（投影到SO(3)流形）
    Rx_ortho = projectToSO3(Rx);
    Ry_ortho = projectToSO3(Ry);
    Rz_ortho = projectToSO3(Rz);
    
    % 组合为齐次变换矩阵
    X = rt2T(Rx_ortho, Tx);
    Y = rt2T(Ry_ortho, Ty);
    Z = rt2T(Rz_ortho, Tz);
end

%%
function error_value = errs(A, B, C, X, Y, Z)
% err - 计算AXB=YCZ校准误差（基于2-范数）
% 输入:
%   A, B, C, X, Y, Z - 4x4齐次变换矩阵
% 输出:
%   error_value - 标量，AXB与YCZ的2-范数差异

    % 输入验证：确保所有输入为4x4矩阵
    if ~isequal(size(A), [4,4]) || ~isequal(size(B), [4,4]) || ~isequal(size(C), [4,4]) || ...
       ~isequal(size(X), [4,4]) || ~isequal(size(Y), [4,4]) || ~isequal(size(Z), [4,4])
        error('所有输入矩阵必须是4x4维度');
    end
    
    % 计算变换链：A*X*B 和 Y*C*Z
    AXB = A * X * B;
    YCZ = Y * C * Z;
    
    % 计算差异矩阵
    difference = AXB - YCZ;
    
    % 计算2-范数（谱范数）
    error_value = norm(difference, 2);
end

%%
function error_val = errfunc_wang(V, P, grad_V, eta, conf)
% errfunc_wang - Wang算法误差函数（基于2-范数）
% 输入:
%   V - 优化变量，通常为3x12矩阵 [Rx,Ry,Rz,Tx,Ty,Tz]
%   P - 元胞数组，包含校准数据 {A_cell, B_cell, C_cell}
%       A_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人1基座到末端法兰）
%       B_cell: 元胞数组，每个元素为4x4齐次变换矩阵（工具或传感器变换）
%       C_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人2末端法兰到工具）
%   grad_V - 梯度向量（未使用，保留用于兼容性）
%   eta - 学习率（未使用，保留用于兼容性）
%   conf - 配置结构体（未使用，保留用于兼容性）
% 输出:
%   error_val - 标量，平均误差值

    % 输入验证
    if ~iscell(P) || length(P) ~= 3
        error('输入P必须是包含3个元胞数组的细胞数组');
    end
    
    % 从优化变量V提取变换矩阵X, Y, Z
    [X, Y, Z] = post_process_Wang(V);
    
    % 提取数据
    A_cell = P{1};
    B_cell = P{2};
    C_cell = P{3};
    
    % 验证数据一致性
    n = length(A_cell);
    if length(B_cell) ~= n || length(C_cell) ~= n
        error('A、B、C数据点数量必须一致');
    end
    
    % 计算所有数据点的总误差
    total_error = 0;
    for i = 1:n
        % 获取当前数据点
        A_i = A_cell{i};
        B_i = B_cell{i};
        C_i = C_cell{i};
        
        % 计算当前数据点的误差
        error_i = errs(A_i, B_i, C_i, X, Y, Z);
        total_error = total_error + error_i;
    end
    
    % 计算平均误差
    error_val = total_error / n;
end

%%
function grad = get_one_grad_wang_v(A, B, C, V, mu)
% get_one_grad_wang_v - 计算Wang目标函数对于单个三元组位姿的梯度
% 
% 输入参数:
%   A, B, C: 4x4齐次变换矩阵，表示测量数据中的位姿
%   V: 3x12优化变量矩阵，结构为 [Rx, Ry, Rz, Tx, Ty, Tz]
%      Rx, Ry, Rz为3x3旋转矩阵块，Tx, Ty, Tz为3x1平移向量
%   mu: 5元素权重向量，对应目标函数中各项的权重系数
%
% 输出参数:
%   grad: 3x12梯度矩阵，包含目标函数对优化变量的偏导数
%
% 算法对应: 
%   Wang et al. "Simultaneous Calibration of Multicoordinates for a Dual-Robot System by Solving the AXB=YCZ Problem"
%   中的梯度计算公式(51)

    % 输入验证
    if size(A,1) ~= 4 || size(A,2) ~= 4 || ...
       size(B,1) ~= 4 || size(B,2) ~= 4 || ...
       size(C,1) ~= 4 || size(C,2) ~= 4
        error('输入A、B、C必须是4x4齐次变换矩阵');
    end
    
    if size(V,1) ~= 3 || size(V,2) ~= 12
        error('输入V必须是3x12优化变量矩阵');
    end
    
    if length(mu) ~= 5
        error('权重向量mu必须包含5个元素');
    end

    % 从优化变量V中提取旋转矩阵和平移向量
    Rx = V(:, 1:3);   % 提取Rx（3x3矩阵）
    Ry = V(:, 4:6);   % 提取Ry（3x3矩阵）
    Rz = V(:, 7:9);   % 提取Rz（3x3矩阵）
    Tx = V(:, 10);    % 提取Tx（3x1向量）
    Ty = V(:, 11);    % 提取Ty（3x1向量）
    Tz = V(:, 12);    % 提取Tz（3x1向量）

    % 从变换矩阵A、B、C中提取旋转部分和平移部分
    Ra = A(1:3, 1:3);  % A的旋转矩阵部分
    Rb = B(1:3, 1:3);  % B的旋转矩阵部分
    Rc = C(1:3, 1:3);  % C的旋转矩阵部分
    Ta = A(1:3, 4);    % A的平移向量部分
    Tb = B(1:3, 4);    % B的平移向量部分
    Tc = C(1:3, 4);    % C的平移向量部分

    % 计算梯度分量1：对应旋转矩阵Rx的梯度
    term1_rx = mu(1) * (Rx * Rb * Rb' - Ra' * Ry * Rc * Rz * Rb');
    term2_rx = mu(2) * (Rx * (Tb * Tb') + Tx * Tb' + Ra' * Ta * Tb' - ...
               Ra' * Ry * Rc * Tz * Tb' - Ra' * Ry * Tc * Tb' - Ra' * Ty * Tb');
    term3_rx = 2 * mu(3) * (Rx * Rx' * Rx - Rx);
    row1 = term1_rx + term2_rx + term3_rx;

    % 计算梯度分量2：对应旋转矩阵Ry的梯度
    term1_ry = mu(1) * (Ry * Rc * Rz * Rz' * Rc' - Ra * Rx * Rb * Rz' * Rc');
    term2_ry = mu(2) * (Ry * (Tc * Tc') + Ry * Rc * (Tz * Tz') * Rc' - ...
               Ra * Rx * Tb * Tz' * Rc' - Ra * Rx * Tb * Tc' - Ra * Tx * Tz' * Rc' - ...
               Ra * Tx * Tc' - Ta * Tz' * Rc' - Ta * Tc' + Ry * Tc * Tz' * Rc' + ...
               Ry * Rc * Tz * Tc' + Ty * Tz' * Rc' + Ty * Tc');
    term3_ry = 2 * mu(4) * (Ry * Ry' * Ry - Ry);
    row2 = term1_ry + term2_ry + term3_ry;

    % 计算梯度分量3：对应旋转矩阵Rz的梯度
    term1_rz = mu(1) * (Rc' * Ry' * Ry * Rc * Rz - Rc' * Ry' * Ra * Rx * Rb);
    term2_rz = 2 * mu(5) * (Rz * Rz' * Rz - Rz);
    row3 = term1_rz + term2_rz;

    % 计算梯度分量4：对应平移向量Tx的梯度
    row4 = mu(2) * (Tx + Rx * Tb + Ra' * Ta - Ra' * Ry * Rc * Tz - Ra' * Ry * Tc - Ra' * Ty);

    % 计算梯度分量5：对应平移向量Ty的梯度
    row5 = mu(2) * (Ty - Ra * Rx * Tb - Ra * Tx - Ta + Ry * Rc * Tz + Ry * Tc);

    % 计算梯度分量6：对应平移向量Tz的梯度
    row6 = mu(2) * (Rc' * Ry' * Ry * Rc * Tz - Rc' * Ry' * Ra * Rx * Tb - ...
           Rc' * Ry' * Ra * Tx - Rc' * Ry' * Ta + Rc' * Ry' * Ry * Tc + Rc' * Ry' * Ty);

    % 组合所有梯度分量形成3x12梯度矩阵
    grad = [row1, row2, row3, row4, row5, row6];
end

%%
function grad = get_i_grad_wang_v(V, P, i, conf)
% get_i_grad_wang_v - 计算第i个数据点的Wang目标函数梯度
% 
% 输入参数:
%   V - 3x12优化变量矩阵，结构为 [Rx, Ry, Rz, Tx, Ty, Tz]
%   P - 元胞数组，包含校准数据 {A_cell, B_cell, C_cell}
%       A_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人1基座到末端法兰）
%       B_cell: 元胞数组，每个元素为4x4齐次变换矩阵（工具或传感器变换）
%       C_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人2末端法兰到工具）
%   i - 标量，数据点索引（从1开始）
%   conf - 结构体，包含配置参数，必须包含mu字段（5元素权重向量）
%
% 输出参数:
%   grad - 3x12梯度矩阵，目标函数对优化变量的偏导数

    % 输入验证
    if ~iscell(P) || length(P) ~= 3
        error('输入P必须是包含3个元胞数组的细胞数组');
    end
    
    A_cell = P{1};
    B_cell = P{2};
    C_cell = P{3};
    
    % 验证索引有效性
    n = length(A_cell);
    if i < 1 || i > n
        error('索引i=%d超出范围，有效范围: 1-%d', i, n);
    end
    
    if length(B_cell) ~= n || length(C_cell) ~= n
        error('A、B、C数据点数量不一致: A=%d, B=%d, C=%d', ...
              length(A_cell), length(B_cell), length(C_cell));
    end
    
    % 验证配置参数
    if ~isfield(conf, 'mu')
        error('配置结构体conf必须包含mu字段');
    end
    
    if length(conf.mu) ~= 5
        error('权重向量mu必须包含5个元素');
    end
    
    % 提取第i个数据点
    A_i = A_cell{i};
    B_i = B_cell{i};
    C_i = C_cell{i};
    
    % 验证数据点格式
    if size(A_i,1) ~= 4 || size(A_i,2) ~= 4 || ...
       size(B_i,1) ~= 4 || size(B_i,2) ~= 4 || ...
       size(C_i,1) ~= 4 || size(C_i,2) ~= 4
        error('数据点A、B、C必须是4x4齐次变换矩阵');
    end
    
    % 调用核心梯度计算函数
    grad = get_one_grad_wang_v(A_i, B_i, C_i, V, conf.mu);
end

%%
function grad_full = get_full_grad_wang_v(V, P, conf)
% get_full_grad_wang_v - 计算AXB=YCZ问题目标函数对整个数据集的平均梯度
%
% 输入参数:
%   V - 3x12优化变量矩阵，结构为 [Rx, Ry, Rz, Tx, Ty, Tz]
%       其中Rx, Ry, Rz为3x3旋转矩阵块，Tx, Ty, Tz为3x1平移向量
%   P - 元胞数组，包含校准数据集 {A_cell, B_cell, C_cell}
%       A_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人1基座到末端法兰）
%       B_cell: 元胞数组，每个元素为4x4齐次变换矩阵（工具或传感器变换）  
%       C_cell: 元胞数组，每个元素为4x4齐次变换矩阵（机器人2末端法兰到工具）
%   conf - 结构体，包含配置参数，必须包含mu字段（5元素权重向量）
%
% 输出参数:
%   grad_full - 3x12平均梯度矩阵，目标函数对优化变量的偏导数平均值
%
% 算法对应:
%   Wang et al. "Simultaneous Calibration of Multicoordinates for a Dual-Robot System"
%   中的全梯度计算公式，用于SVRG等优化算法

    % 输入验证
    if ~iscell(P) || length(P) ~= 3
        error('输入P必须是包含3个元胞数组的细胞数组 {A_data, B_data, C_data}');
    end
    
    % 从P中提取数据集合
    A_cell = P{1};
    B_cell = P{2};
    C_cell = P{3};
    
    % 获取数据点数量并验证一致性
    n = length(A_cell);
    if length(B_cell) ~= n || length(C_cell) ~= n
        error('A、B、C数据集大小不一致: A=%d, B=%d, C=%d', n, length(B_cell), length(C_cell));
    end
    
    if n == 0
        error('数据集不能为空');
    end
    
    % 验证配置参数
    if ~isfield(conf, 'mu') || length(conf.mu) ~= 5
        error('配置必须包含5元素权重向量conf.mu');
    end
    
    % 初始化梯度累加器
    grad_full = zeros(size(V));  % 3x12零矩阵
    
    % 循环计算每个数据点的梯度并累加
    for i = 1:n
        % 获取当前数据点
        A_i = A_cell{i};
        B_i = B_cell{i};
        C_i = C_cell{i};
        
        % 验证数据点格式
        if size(A_i,1) ~= 4 || size(A_i,2) ~= 4 || ...
           size(B_i,1) ~= 4 || size(B_i,2) ~= 4 || ...
           size(C_i,1) ~= 4 || size(C_i,2) ~= 4
            error('数据点%d的变换矩阵必须是4x4齐次变换矩阵', i);
        end
        
        % 计算当前数据点的梯度
        grad_i = get_one_grad_wang_v(A_i, B_i, C_i, V, conf.mu);
        
        % 累加到总梯度
        grad_full = grad_full + grad_i;
    end
    
    % 计算平均梯度
    grad_full = grad_full / n;
end

%%
function x_new = update_point(y_t, y_t1, x_t1, tau)
% update_point - Katyusha算法的点更新规则
% 输入:
%   y_t, y_t1 - 当前和前一次迭代点
%   x_t1 - 前一次快照点
%   tau - 动量参数
% 输出:
%   x_new - 更新后的快照点

    % Katyusha特定的线性组合更新
    % 对应Julia代码中的更新公式: (3/2*y_t + 1/2*x_t1 - (1-tau)*y_t1) / (1+tau)
    numerator = (3.0/2)*y_t + (1.0/2)*x_t1 - (1 - tau)*y_t1;
    x_new = numerator / (1.0 + tau);
end

%%
function [y_t, err_history] = katyushaX(V, P, get_i_grad_func, get_full_grad_func, err_func, normfunc, conf)
% katyushaX - Katyusha优化算法实现，用于AXB=YCZ校准问题
%
% 输入参数:
%   V - 3x12优化变量矩阵，初始估计值
%   P - 元胞数组，包含校准数据 {A_cell, B_cell, C_cell}
%   get_i_grad_func - 函数句柄，计算单个数据点的梯度
%   get_full_grad_func - 函数句柄，计算全数据集的梯度
%   err_func - 函数句柄，计算误差
%   normfunc - 函数句柄，对优化变量进行归一化/投影
%   conf - 结构体，包含算法配置参数
%
% 输出参数:
%   y_t - 3x12矩阵，优化后的变量
%   err_history - 误差历史记录，用于收敛分析

    % 从配置中提取参数
    max_iter = conf.max_iter;      % 最大迭代次数
    m = conf.m;                    % 内部迭代次数
    eta0 = conf.eta;               % 学习率
    n = length(P{1});              % 数据点数量
    tau = conf.tau;                % 动量参数
    
    % 定义梯度计算函数
    grad_full = @(V) get_full_grad_func(V, P, conf);
    grad_i = @(V, i) get_i_grad_func(V, P, i, conf);
    
    % 初始化变量
    y_t = V;          % 当前迭代点
    y_t1 = V;         % 前一次迭代点
    x_t = V;          % 快照点
    x_t1 = V;         % 前一次快照点
    
    % 计算初始全梯度
    grad_V_tilde = grad_full(x_t);
    
    % 初始化误差记录
    if isfield(conf, 'err') && strcmp(conf.err, 'Ext')
        err_history = zeros(max_iter, conf.errdim);
    else
        err_history = zeros(max_iter, 1);
    end
    
    % 主优化循环
    for iter = 1:max_iter
        % 保存前一次快照点
        x_t1 = x_t;
        
        % 更新快照点（Katyusha特定更新规则）
        x_t = update_point(y_t, y_t1, x_t1, tau);
        
        % 计算快照点的全梯度（用于方差缩减）
        grad_V_tilde = grad_full(x_t);
        
        % 内部迭代（随机方差缩减）
        y_k = x_t;  % 从快照点开始内部迭代
        
        for k = 1:m
            % 随机选择一个数据点
            i_k = randi(n);
            
            % Katyusha更新规则：结合当前点梯度和快照点梯度
            grad_current = grad_i(y_k, i_k);
            grad_snapshot = grad_i(x_t, i_k);
            
            % 梯度更新
            gradient_update = grad_current - grad_snapshot + grad_V_tilde;
            y_k = y_k - eta0 * gradient_update;
        end
        
        % 更新迭代点
        y_t1 = y_t;
        
        % 应用正则化或归一化
        if isfield(conf, 'reg') && conf.reg
            y_t = y_k;  % 不使用归一化
        else
            y_t = normfunc(y_k);  % 应用归一化/投影
        end
        
        % 记录误差
        if isfield(conf, 'err') && strcmp(conf.err, 'Ext')
            err_history(iter, :) = err_func(y_t, P, grad_V_tilde, 0.0, conf);
        else
            err_history(iter) = err_func(y_t, P, grad_V_tilde, 0.0, conf);
        end
        
        % 收敛检查（可选）
        if isfield(conf, 'tolerance') && err_history(iter) < conf.tolerance
            fprintf('收敛于迭代 %d，误差: %.6f\n', iter, err_history(iter));
            err_history = err_history(1:iter);  % 截断误差历史
            break;
        end
        
        % 进度显示
        if isfield(conf, 'verbose') && conf.verbose && mod(iter, 100) == 0
            fprintf('迭代 %d: 误差 = %.6f\n', iter, err_history(iter));
        end
    end
end


%% 
function [Result, Err] = iter_solve_Katyusha(V, P, conf)
% iter_solve_Katyusha - 使用Katyusha优化算法进行迭代优化
% 输入:
%   V - 元胞数组，包含3个4x4齐次变换矩阵的初始估计 {X_init, Y_init, Z_init}
%   P - 元胞数组，包含校准数据 {A_cell, B_cell, C_cell}，每个是4x4变换矩阵的元胞数组
%   conf - 结构体，包含配置参数
% 输出:
%   Result - 优化后的变量（3x12矩阵或优化后的变换矩阵）
%   Err - 误差历史记录

    % 将初始变换矩阵转换为优化变量格式
    V_init = toV(V);
    
    % 根据配置选择误差函数
    if isfield(conf, 'err')
        switch conf.err
            case 'NAN'
                errfunc = @errN;  % 空误差函数
            case 'LS'
                errfunc = @errfunc_wang;  % 最小二乘误差
            case 'Ext'
                if isfield(conf, 'errfunc')
                    errfunc = conf.errfunc;  % 外部提供的误差函数
                else
                    error('外部误差函数未提供');
                end
            otherwise
                error('不支持的误差类型: %s', conf.err);
        end
    else
        % 默认使用最小二乘误差
        errfunc = @errfunc_wang;
    end
    
    % 调用Katyusha优化算法
    [Result, Err] = katyushaX(V_init, P, @get_i_grad_wang_v, @get_full_grad_wang_v, errfunc, @normfunc, conf);
end
