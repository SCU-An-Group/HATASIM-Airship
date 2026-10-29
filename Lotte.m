% 【Lotte】Last modified: 2026.6.26
clc; clear; close all;    % 清空命令行、工作区、图表

%% ======== ★参数设置★ ========
% 浮空器参数
V_ast = 109;                  % 体积 (m³)
L_ast = V_ast^(1/3);          % 参考长度 (m)
S_ref = V_ast^(2/3);          % 参考面积 (m²)【体积当量面积】
m_enve = 100;                 % 囊体质量 (kg)
m_load = 10;                  % 载荷质量 (kg)
x_CB = 0;    y_CB = 0;        % 浮心坐标 (m) 【体轴系：以质心（重心）为原点，X指向机头为正，Y垂直X轴向上为正】

% 系留绳参数
L_seg = 1;                    % 绳段原长 (m)【该值影响迭代均匀性和计算结果精度】
d_tet = 0.004;                % 绳直径 (m)
A_tet = pi/4 * d_tet^2;       % 横截面积 (m²)
rho_tet = 970;                % 绳密度 (kg/m³)
E = 172e9;                    % 弹性模量 (Pa)
x_CP = 0;    y_CP = -4;       % 主系点坐标 (m)【体轴系：以质心（重心）为原点，X指向机头为正，Y垂直X轴向上为正】

% 风场参数（幂律模型）
v_ref = 2;                    % 参考风速 (m/s)【设为0即无风】
h_ref = 10;                   % 参考风速高度 (m)
h_ground = 0;                 % 地面海拔高度 (m)
alpha_ws = 0.143;             % 风切变指数【科普中国默认1/7、设为0即定常风】

% 大气与物理常量
R_earth = 6356766;            % 地球半径 (m)
g0 = 9.80665;                 % 海平面重力加速度 (m/s²)
P0 = 101325;                  % 海平面气压 (Pa)
T0 = 288.15;                  % 海平面温度 (K)
gamma = 0.0065;               % 温度垂直递减率 (K/m)
R_air = 287.05287;            % 空气气体常数 (J/(kg·K))
R_He = 2077.2298;             % 氦气气体常数 (J/(kg·K))
mu0 = 1.7894e-5;              % 海平面动力粘度 (Pa·s)
S = 110.4;                    % 萨瑟兰常数 (K)

%% ======== 初始状态多解遍历 ========
% 计算第一个有效高度及其参数
h0 = h_ground + abs(y_CP);                    % 系留点刚离开地面时的海拔高度
[g0_h, T0_h, P0_h, rho_air0, rho_He0, mu0_h, v0] = ISA_params( ...
    h0, R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, v_ref, h_ref, alpha_ws);
m_He0 = rho_He0 * V_ast;                      % 初始氦气质量
W_ast0 = (m_enve + m_load + m_He0) * g0_h;    % 初始浮空器总重
F_b0 = rho_air0 * V_ast * g0_h;               % 初始浮力

% 存放有效解
valid_solutions = [];     % 每行 [alpha, T_CP, theta]
valid_init_alpha = [];    % 对应的初始猜测α
fprintf('======== 正在求解初始高度 h = %.2f m 的全部有效解 ========\n', h0);
alpha_samples_deg = linspace(-90, 90, 181);        % 度数采样
alpha_samples_rad = deg2rad(alpha_samples_deg);    % 转换为弧度

% 初始迎角遍历求解
for i = 1:length(alpha_samples_rad)
    alpha0 = alpha_samples_rad(i);
    % 调用定高求解
    [alpha, T_CP, theta, ~, ~, valid, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        SolveAtHeight(h0, L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
            d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
            v_ref, h_ref, h_ground, alpha_ws, ...
            R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, alpha0);
    % 判断解是否有效
    if valid
        if alpha >= -pi/6 && alpha <= pi/6 && T_CP >= 0 && theta >= 0 && theta <= pi/2    % ★有效解范围★
            valid_solutions = [valid_solutions; alpha, T_CP, theta];
            valid_init_alpha = [valid_init_alpha; alpha0];
        end
    end
end

% α解去重
tol_alpha = deg2rad(0.01);      % 设置α解的容差 (°)
if ~isempty(valid_solutions)    % 有解存储
    [~, ia] = uniquetol(valid_solutions(:,1), tol_alpha, 'DataScale', 1);
    unique_solutions = valid_solutions(ia, :);
    unique_init_alpha = valid_init_alpha(ia, :);
else
    unique_solutions = [];    unique_init_alpha = [];
end
if isempty(unique_solutions)    % 无解报错
    error('在初始高度 h = %.2f m 未找到任何有效解，请检查参数！', h0);
end

% 输出所有独立有效解
fprintf('\n共找到【%d】组独立有效解：\n', size(unique_solutions,1));
for i = 1:size(unique_solutions,1)
    fprintf('解 %d: 初始迎角 = %.2f°，初始拉力 = %.2f N，初始拉力夹角 = %.2f°\n', ...
            i, unique_solutions(i,1)*180/pi, unique_solutions(i,2), unique_solutions(i,3)*180/pi);
end

%% ======== 自适应步进迭代函数 ========
function [H_results, X_results, F_b_results, W_ast_results, Alpha_results, ...
          F_l_results, F_d_results, M_CG_results, T_CP_results, Theta_results, ...
          L_tet_ori_results, L_tet_def_results, W_tet_results, final_tet_data, ...
          final_z_CP, final_x_CP, final_F_ground, final_g, final_T, final_P, ...
          final_rho_air, final_rho_He, final_mu, final_v, init_solution_at_h0, ...
          dM_dv_results, dM_da_results, stress_all_min, stress_all_max] = ...
    ComputeTrajectory(start_h, init_alpha, ...
                      L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
                      d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
                      v_ref, h_ref, h_ground, alpha_ws, ...
                      R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S)
    % 初始化结果数组
    H_results = [];            X_results = [];            F_b_results = [];    W_ast_results = [];
    Alpha_results = [];        F_l_results = [];          F_d_results = [];    M_CG_results = [];
    dM_dv_results = [];        dM_da_results = [];        T_CP_results = [];   Theta_results = [];
    L_tet_ori_results = [];    L_tet_def_results = [];    W_tet_results = [];        
    stress_all_min = inf;      stress_all_max = -inf;
    
    % 迭代求解参数
    max_iter = 10000;                       % 最大迭代次数
    h_step_init = V_ast^(1/3);              % 初始步进高度 (m)
    h_step = h_step_init;                   % 自适应步进高度 (m)
    delta_S = V_ast^(1/3);                  % 升空轨迹迭代间距 (m)
    current_h = start_h - h_step;           % 从起点下方一个步长处开始【使得第一个点刚好是初始点】
    last_alpha = init_alpha;                % 起始高度的初始α

    % 迭代求解
    fprintf('======== 正在进行迭代求解 ========\n');
    for iter = 1:max_iter
        if h_step < 0.001
            break;
        end
        attempted_h = current_h + h_step;     % 更新迭代尝试高度
        
        % 调用定高求解（传入上次α解作为初值）
        [alpha, T_CP, theta, drift, tet_data, valid, ...
         F_b, W_ast, F_l, F_d, M_CG, L_tet_ori, L_tet_def, W_tet, ...
         g, T, P, rho_air, rho_He, mu, v, z_CP, x_cum, F_ground, ...
         dM_dv, dM_da] = SolveAtHeight(...
            attempted_h, L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
            d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
            v_ref, h_ref, h_ground, alpha_ws, ...
            R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, last_alpha);

        if valid
            % 记录随高度变化的数据
            H_results = [H_results, attempted_h];
            X_results = [X_results, drift];
            F_b_results = [F_b_results, F_b];
            W_ast_results = [W_ast_results, W_ast];
            Alpha_results = [Alpha_results, alpha];
            F_l_results = [F_l_results, F_l];
            F_d_results = [F_d_results, F_d];
            M_CG_results = [M_CG_results, M_CG];
            T_CP_results = [T_CP_results, T_CP];
            Theta_results = [Theta_results, theta];
            L_tet_ori_results = [L_tet_ori_results, L_tet_ori];
            L_tet_def_results = [L_tet_def_results, L_tet_def];
            W_tet_results = [W_tet_results, W_tet];
            dM_dv_results = [dM_dv_results, dM_dv];
            dM_da_results = [dM_da_results, dM_da];

            % 更新全局应力极值
            if ~isempty(tet_data)
                stress_h = tet_data(:,5);
                stress_h_min = min(stress_h);
                stress_h_max = max(stress_h);
                if stress_h_min < stress_all_min; stress_all_min = stress_h_min; end
                if stress_h_max > stress_all_max; stress_all_max = stress_h_max; end
            end

            % 更新最终高度下的数据
            final_tet_data = tet_data;
            final_z_CP = z_CP;
            final_x_CP = x_cum;
            final_F_ground = F_ground;
            final_g = g;
            final_T = T;
            final_P = P;
            final_rho_air = rho_air;
            final_rho_He = rho_He;
            final_mu = mu;
            final_v = v;
            current_h = attempted_h;    % 更新当前高度
            last_alpha = alpha;         % 记录本次结果作为下次初始值

            % 自适应调整步长（根据前两点的实际弧长）
            if length(H_results) >= 2
                dx = diff(X_results(end-1:end));
                dh = diff(H_results(end-1:end));
                ds = sqrt(dx^2 + dh^2);
                dh_ds = abs(dh / (ds + 1e-6));
                h_step = delta_S * dh_ds;
            else
                h_step = h_step_init;
            end
        else
            h_step = h_step * 0.5;    % 末端步幅缩进提升结果精度
        end
    end

    % 检查是否存在有效解
    if isempty(H_results)
        error('未找到任何有效高度，请检查参数设置！');
    end
    
    % 从第一个结果中提取起点解
    init_solution_at_h0 = [Alpha_results(1), T_CP_results(1), Theta_results(1)];
end

%% 选择解进行迭代与交互
while true
    fprintf('\n=========================================================\n');
    prompt = sprintf('请选择要运行的解编号 (1~%d)，输入 0 或直接回车退出程序：\n请输入编号：', size(unique_solutions,1));
    choice = input(prompt);
    if isempty(choice) || choice == 0
        break;
    end
    if ~isscalar(choice) || choice < 1 || choice > size(unique_solutions,1)
        fprintf('输入无效，请重新输入！\n');
        continue;
    end

    close all;    % 关闭之前的所有图形窗口，重新绘制
    init_alpha_h0 = unique_solutions(choice, 1);    % 获取选中解的初始α
    start_h = h0;    % 起始高度

    % 调用步进迭代函数
    [H_results, X_results, F_b_results, W_ast_results, Alpha_results, ...
     F_l_results, F_d_results, M_CG_results, T_CP_results, Theta_results, ...
     L_tet_ori_results, L_tet_def_results, W_tet_results, final_tet_data, ...
     final_z_CP, final_x_CP, final_F_ground, final_g, final_T, final_P, ...
     final_rho_air, final_rho_He, final_mu, final_v, init_solution_at_h0, ...
     dM_dv_results, dM_da_results, stress_all_min, stress_all_max] = ...
        ComputeTrajectory(start_h, init_alpha_h0, ...
                          L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
                          d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
                          v_ref, h_ref, h_ground, alpha_ws, ...
                          R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S);

    % 对结果排序（按高度升序）
    [H_results_sorted, sort_idx] = sort(H_results);
    Alpha_sorted = Alpha_results(sort_idx) * 180/pi;
    T_CP_sorted = T_CP_results(sort_idx);

    % 构建插值函数
    alpha_interp = griddedInterpolant(H_results_sorted, Alpha_sorted, 'pchip', 'linear');
    T_interp = griddedInterpolant(H_results_sorted, T_CP_sorted, 'pchip', 'linear');
    h_min = H_results_sorted(1);    h_max = H_results_sorted(end);    % 最小/大高度

    % 求最大迎角
    opt_alpha = @(h) -alpha_interp(h);
    [h_alpha_max, neg_alpha_max] = fminbnd(opt_alpha, h_min, h_max);
    alpha_max = -neg_alpha_max;

    % 求最大拉力
    opt_T = @(h) -T_interp(h);
    [h_T_max, neg_T_max] = fminbnd(opt_T, h_min, h_max);
    T_max = -neg_T_max;
    stress_max_from_T = (T_max / A_tet) * (1 + T_max / (E * A_tet)) / 1e6;   % MPa

    % 输出最终结果
    fprintf('\n--------【解 %d 浮空器最终结果】--------\n', choice);
    fprintf('最大驻空高度: %.2f m\n', H_results(end));
    fprintf('最远漂移距离: %.2f m\n', X_results(end));
    fprintf('迎角: %.2f° (最大值: %.2f° @ %.2f m)\n', Alpha_results(end)*180/pi, alpha_max, h_alpha_max);
    fprintf('气动升力: %.2f N\n', F_l_results(end));
    fprintf('气动阻力: %.2f N\n', F_d_results(end));
    fprintf('俯仰力矩: %.2f N·m\n', M_CG_results(end));
    fprintf('浮空器浮力: %.2f N\n', F_b_results(end));
    fprintf('浮空器重力: %.2f N\n', W_ast_results(end));
    fprintf('动压静稳定性导数: %.4f N·m/(m/s)\n', dM_dv_results(end));
    fprintf('迎角静稳定性导数: %.2f N·m/rad (%.2f N·m/°)\n', dM_da_results(end), dM_da_results(end)*pi/180);

    fprintf('\n--------【解 %d 系留绳最终结果】--------\n', choice);
    fprintf('系留绳原长: %.2f m\n', L_tet_ori_results(end));
    fprintf('系留绳形变后长度: %.2f m\n', L_tet_def_results(end));
    fprintf('主系点拉力: %.2f N (最大值: %.2f N, %.2f MPa @ %.2f m)\n', ...
    T_CP_results(end), T_max, stress_max_from_T, h_T_max);
    fprintf('主系点拉力夹角: %.2f°\n', Theta_results(end)*180/pi);
    fprintf('系留绳重力: %.2f N\n', W_tet_results(end));

    fprintf('\n--------【解 %d 大气环境最终参数】--------\n', choice);
    fprintf('风速: %.2f m/s\n', final_v);
    fprintf('温度: %.2f K (%.2f °C)\n', final_T, final_T-273.15);
    fprintf('压强: %.2f Pa\n', final_P);
    fprintf('空气密度: %.4f kg/m³\n', final_rho_air);
    fprintf('氦气密度: %.4f kg/m³\n', final_rho_He);
    fprintf('动力粘度: %.4e Pa·s\n', final_mu);
    fprintf('重力加速度: %.4f m/s²\n', final_g);

    fprintf('\n--------【解 %d 静稳定性分析】--------\n', choice);
    fprintf('动压静稳定性导数: %.4f N·m/(m/s)\n', dM_dv_results(end));
    fprintf('迎角静稳定性导数: %.2f N·m/rad (%.2f N·m/°)\n', dM_da_results(end), dM_da_results(end)*pi/180);
    if all(dM_dv_results < 0)    % 动压静稳定性导数小于0则稳定
        fprintf('是否全过程具有动压静稳定性（迎风低头）：是\n');
    else
        fprintf('是否全过程具有动压静稳定性（迎风低头）：否\n');
    end
    if all(dM_da_results < 0)    % 迎角静稳定性导数小于0则稳定
        fprintf('是否全过程具有迎角静稳定性（自动回正）：是\n');
    else
        fprintf('是否全过程具有迎角静稳定性（自动回正）：否\n');
    end

    fprintf('\n======== 迭代求解完毕！========\n\n');

    %% ======== 结果可视化 ========
    % -------- 图一：升空轨迹与系留绳形态 --------
    figure('Name', '升空轨迹与系留绳形态');

    % 绘制升空轨迹和最终点
    traj_handle = plot(X_results, H_results, 'k-', 'LineWidth', 3, 'DisplayName', '升空轨迹');
    trail_after_handle = [];   % 初始化查询点之后的点线轨迹句柄
    hold on;
    plot(X_results(end), H_results(end), 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', sprintf('最终点'));

    % 提取系留绳数据
    n_seg = size(final_tet_data, 1);     % 绳段数量
    seg_F = final_tet_data(:,1);         % 绳段张力
    seg_theta = final_tet_data(:,2);     % 张力夹角
    seg_z_down = final_tet_data(:,3);    % 绳段下端点高度
    seg_L_def = final_tet_data(:,4);     % 绳段形变后长度
    seg_stress = final_tet_data(:,5);    % 绳段应力

    % 计算微元段端点坐标
    x_nodes = zeros(n_seg+1, 1);
    z_nodes = zeros(n_seg+1, 1);
    x_nodes(1) = final_x_CP;                      % 主系点水平坐标
    z_nodes(1) = final_z_CP;                      % 主系点海拔高度
    for i = 1:n_seg
        dx = seg_L_def(i) * sin(seg_theta(i));    % 该段水平投影长度
        x_nodes(i+1) = x_nodes(i) - dx;           % 下一节点水平坐标
        z_nodes(i+1) = seg_z_down(i);             % 下一段上端点高度 = 上一段下端点高度
    end
    x_nodes(end) = 0;                             % 强制锚点横坐标为0，消除数值累积误差

    % 绘制系留绳形态（按应力着色），使用全局应力范围
    if stress_all_max > stress_all_min
        cmap = colormap(turbo(256));                   
        stress_norm = (seg_stress - stress_all_min) / (stress_all_max - stress_all_min);
        idx = max(1, min(size(cmap,1), round(stress_norm * (size(cmap,1)-1) + 1)));
        for i = 1:n_seg
            seg_x = [x_nodes(i), x_nodes(i+1)];
            seg_y = [z_nodes(i), z_nodes(i+1)];
            line(seg_x, seg_y, 'Color', cmap(idx(i),:), 'LineWidth', 4, 'HandleVisibility', 'off');
        end
    else
        % 所有应力相同，使用单一颜色
        for i = 1:n_seg
            seg_x = [x_nodes(i), x_nodes(i+1)];
            seg_y = [z_nodes(i), z_nodes(i+1)];
            line(seg_x, seg_y, 'Color', [0 0 1], 'LineWidth', 4, 'HandleVisibility', 'off');
        end
    end
    
    % 设置图一格式
    xlabel('Drift, X [m]');
    ylabel('Altitude, H [m]');
    title('升空轨迹与系留绳形态');
    legend('Location', 'best');
    axis equal;
    grid on;
    c = colorbar;
    c.Label.String = 'Stress, [MPa]';
    clim([stress_all_min, stress_all_max]);

    % -------- 图二：系留绳张力分布 --------
    figure('Name', '系留绳张力分布');

    % 计算每个绳段参数
    L_def_seg = seg_L_def;                                     % 形变后长度
    cum_L_from_CP = cumsum([0; L_def_seg]);                    % 距主系点的长度
    L_from_ground = L_tet_def_results(end) - cum_L_from_CP;    % 距地面的长度
    T_nodes = [seg_F; final_F_ground];                         % 节点张力

    % 绘制主系点拉力
    T_CP_handle = plot(L_tet_def_results, T_CP_results, 'k-', 'LineWidth', 3, 'DisplayName', '主系点拉力');
    T_CP_after_handle = [];    % 初始化查询点之后的点线句柄
    hold on;
    plot(L_from_ground, T_nodes, 'r-', 'LineWidth', 3, 'DisplayName', '最终绳内张力');

    % 设置图二格式
    xlabel('Length From the Ground, [m]');
    ylabel('Tension, T [N]');
    title('系留绳张力分布');
    legend('Location', 'best');
    grid on;

    %% ======== 参数关系绘制 ========
    % 定义所有可绘制参数的数据、名称及单位
    params_results = {H_results, X_results, F_b_results, W_ast_results, Alpha_results*180/pi, ...
                      F_l_results, F_d_results, M_CG_results, T_CP_results, Theta_results*180/pi, ...
                      L_tet_ori_results, L_tet_def_results, W_tet_results, ...
                      dM_dv_results, dM_da_results};
    params_names = {'浮空器高度 (m)', '浮空器漂移距离 (m)', '浮空器浮力 (N)', '浮空器重力 (N)', '浮空器迎角 (°)', ...
                    '浮空器升力 (N)', '浮空器阻力 (N)', '浮空器俯仰力矩 (N·m)', '主系点拉力 (N)', ...
                    '主系点拉力夹角 (°)', '系留绳原长 (m)', '系留绳形变后长度 (m)', '系留绳重力 (N)', ...
                    '动压静稳定性导数 (N·m/(m/s))', '迎角静稳定性导数 (N·m/rad)'};

    % 交互式选择参数
    fprintf('\n========【解 %d 参数关系绘制】========\n', choice);
    for i = 1:length(params_names)
        fprintf('  %2d : %s\n', i, params_names{i});
    end

    while true
        x_idx = input('\n请输入X轴参数编号 (1-15)，输入0或直接回车结束: ');
        if isempty(x_idx) || x_idx == 0
            break;
        end
        if ~isscalar(x_idx) || x_idx < 1 || x_idx > 15
            fprintf('输入无效，请重新输入1~15之间的整数！\n');
            continue;
        end

        y_idx = input('请输入Y轴参数编号 (1-15)，输入0或直接回车结束: ');
        if isempty(y_idx) || y_idx == 0
            break;
        end
        if ~isscalar(y_idx) || y_idx < 1 || y_idx > 15
            fprintf('输入无效，请重新输入1~15之间的整数！\n');
            continue;
        end

        % 提取绘图数据
        x_data = params_results{x_idx};
        y_data = params_results{y_idx};
        if length(x_data) ~= length(y_data)
            fprintf('数据长度不一致，无法绘图！\n');
            continue;
        end

        % 按迭代顺序绘图
        figure('Name', sprintf('参数关系：%s & %s', params_names{x_idx}, params_names{y_idx}));
        plot(x_data, y_data, 'b-', 'LineWidth', 3);
        hold on;
        plot(x_data(1), y_data(1), 'go', 'MarkerSize', 9, 'MarkerFaceColor', 'g');        % 标记起点
        plot(x_data(end), y_data(end), 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');    % 标记终点

        % 设置参数关系图格式
        xlabel(params_names{x_idx});
        ylabel(params_names{y_idx});
        title(sprintf('%s 随 %s 变化曲线', params_names{y_idx}, params_names{x_idx}));
        legend('变化曲线', '起点', '终点', 'Location', 'best');
        grid on;
        hold off;
        fprintf('图形已绘制完成。\n');
    end
    fprintf('======== 参数关系绘制结束！========\n\n');

    %% ======== 定高查询 ========
    fprintf('\n========【解 %d 定高查询】========\n', choice);
    % 获取已有图形窗口句柄
    fig1 = findobj('Name', '升空轨迹与系留绳形态');
    fig2 = findobj('Name', '系留绳张力分布');
    if isempty(fig1) || isempty(fig2)
        error('未找到图形窗口，请确保先运行主程序生成图表！');
    end

    % 获取图一的颜色映射和应力颜色范围（此时已使用全局范围）
    figure(fig1);
    cmap = colormap;         % 获取当前颜色映射
    stress_limits = clim;    % 获取当前应力颜色范围

    % 初始化绘图对象句柄
    ast_query = [];    tet_query = [];    T_query = [];    % 查询点浮空器 / 系留绳各段曲线 / 查询点张力分布曲线
    last_query_alpha = init_solution_at_h0(1);             % 获取初始猜测值用于查询

    while true
        fprintf('\n有效高度范围: [%.2f, %.2f] m', h0, max(H_results));
        h_query = input('\n请输入查询高度 (m) [输入0或直接回车退出]: ');    % 输入查询高度
        if isempty(h_query) || h_query == 0
            break;
        end

        % 通过插值获取该高度的迎角估计值作为初值，避免收敛到其他分支
        alpha_guess_deg = alpha_interp(h_query);
        if isnan(alpha_guess_deg) || ~isfinite(alpha_guess_deg)
            alpha_guess = last_query_alpha;    % 插值失败则使用上次值
        else
            alpha_guess = deg2rad(alpha_guess_deg);
        end

        % 调用定高求解函数
        [alpha, T_CP, theta, drift, tet_data, valid, ...
         F_b, W_ast, F_l, F_d, M_CG, L_tet_ori, L_tet_def, W_tet, ...
         g, T, P, rho_air, rho_He, mu, v, z_CP, x_cum, F_ground, ...
         dM_dv, dM_da] = SolveAtHeight(...
            h_query, L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
            d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
            v_ref, h_ref, h_ground, alpha_ws, ...
            R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, alpha_guess);
        if ~valid
            fprintf('高度 %.2f m 无效：浮力不足或系留绳松弛，请重新输入！\n', h_query);
            continue;
        end
        last_query_alpha = alpha;    % 更新初始α，确保查询结果使用同一分支解

        % 输出查询结果
        fprintf('\n--------【解 %d 浮空器 @ %.2f m】--------\n', choice, h_query);
        fprintf('漂移距离: %.2f m\n', drift);
        fprintf('迎角: %.2f°\n', alpha*180/pi);
        fprintf('气动升力: %.2f N\n', F_l);
        fprintf('气动阻力: %.2f N\n', F_d);
        fprintf('俯仰力矩: %.2f N·m\n', M_CG);
        fprintf('浮空器浮力: %.2f N\n', F_b);
        fprintf('浮空器重力: %.2f N\n', W_ast);

        fprintf('\n--------【解 %d 系留绳 @ %.2f m】--------\n', choice, h_query);
        fprintf('系留绳原长: %.2f m\n', L_tet_ori);
        fprintf('系留绳形变后长度: %.2f m\n', L_tet_def);
        fprintf('主系点拉力: %.2f N\n', T_CP);
        fprintf('主系点拉力夹角: %.2f°\n', theta*180/pi);
        fprintf('系留绳重力: %.2f N\n', W_tet);

        fprintf('\n--------【解 %d 大气环境 @ %.2f m】--------\n', choice, h_query);
        fprintf('风速: %.2f m/s\n', v);
        fprintf('温度: %.2f K (%.2f °C)\n', T, T-273.15);
        fprintf('压强: %.2f Pa\n', P);
        fprintf('空气密度: %.4f kg/m³\n', rho_air);
        fprintf('氦气密度: %.4f kg/m³\n', rho_He);
        fprintf('动力粘度: %.4e Pa·s\n', mu);
        fprintf('重力加速度: %.4f m/s²\n', g);
        fprintf('\n================================\n');

        fprintf('\n--------【解 %d 静稳定性分析 @ %.2f m】--------\n', choice, h_query);
        fprintf('动压静稳定性导数: %.4f N·m/(m/s)\n', dM_dv);
        fprintf('迎角静稳定性导数: %.2f N·m/rad (%.2f N·m/°)\n', dM_da, dM_da*pi/180);
        if dM_dv < 0
            fprintf('查询点是否具有动压静稳定性（迎风低头）：是\n');
        else
            fprintf('查询点是否具有动压静稳定性（迎风低头）：否\n');
        end
        if dM_da < 0
            fprintf('查询点是否具有迎角静稳定性（自动回正）：是\n');
        else
            fprintf('查询点是否具有迎角静稳定性（自动回正）：否\n');
        end

        % -------- 更新图一：升空轨迹与系留绳形态 --------
        figure(fig1);

        % 清除上一次查询结果（系留绳和查询点）
        if ~isempty(tet_query)
            for i = 1:length(tet_query)
                if ishandle(tet_query(i))
                    delete(tet_query(i));
                end
            end
        end
        if ~isempty(ast_query) && ishandle(ast_query)
            delete(ast_query);
        end

        % 提取查询点绳段数据
        n_seg = size(tet_data, 1);
        seg_theta = tet_data(:,2);      % 张力夹角
        seg_L_def = tet_data(:,4);      % 形变后长度
        seg_z_down = tet_data(:,3);     % 下端点高度
        seg_stress = tet_data(:,5);     % 应力

        % 计算微元段端点坐标
        x_nodes = zeros(n_seg+1, 1);
        z_nodes = zeros(n_seg+1, 1);
        x_nodes(1) = x_cum;             % 主系点水平坐标
        z_nodes(1) = z_CP;              % 主系点海拔高度
        for i = 1:n_seg
            dx = seg_L_def(i) * sin(seg_theta(i));
            x_nodes(i+1) = x_nodes(i) - dx;
            z_nodes(i+1) = seg_z_down(i);
        end
        x_nodes(end) = 0;               % 强制锚点横坐标为0，消除数值累积误差

        % 按应力着色绘制查询绳段，使用与图一相同的全局应力范围
        if stress_limits(2) > stress_limits(1)
            stress_norm = (seg_stress - stress_limits(1)) / (stress_limits(2) - stress_limits(1));
            idx = max(1, min(size(cmap,1), round(stress_norm * (size(cmap,1)-1) + 1)));
            tet_query = gobjects(n_seg, 1);
            for i = 1:n_seg
                seg_x = [x_nodes(i), x_nodes(i+1)];
                seg_y = [z_nodes(i), z_nodes(i+1)];
                tet_query(i) = line(seg_x, seg_y, 'Color', cmap(idx(i),:), ...
                    'LineWidth', 4, 'HandleVisibility', 'off');
            end
        else
            tet_query = gobjects(n_seg, 1);
            for i = 1:n_seg
                seg_x = [x_nodes(i), x_nodes(i+1)];
                seg_y = [z_nodes(i), z_nodes(i+1)];
                tet_query(i) = line(seg_x, seg_y, 'Color', [0 0 1], ...
                    'LineWidth', 4, 'HandleVisibility', 'off');
            end
        end

        % 绘制查询点位置
        ast_query = plot(drift, h_query, 'bo', 'MarkerSize', 7.5, 'MarkerFaceColor', 'b', ...
                                         'LineWidth', 2, 'DisplayName', '查询点');

        % ====== 更新图一：升空轨迹线型 ======
        idx = find(H_results <= h_query, 1, 'last');
        if ~isempty(idx)
            if idx >= 1
                set(traj_handle, 'XData', X_results(1:idx), 'YData', H_results(1:idx));
                if ~isempty(trail_after_handle) && ishandle(trail_after_handle)
                    delete(trail_after_handle);
                end
                if idx < length(H_results)
                    trail_after_handle = plot(X_results(idx+1:end), H_results(idx+1:end), ...
                        'k:', 'LineWidth', 3, 'DisplayName', '潜在升空轨迹');
                else
                    trail_after_handle = [];
                end
            else
                set(traj_handle, 'XData', [], 'YData', []);
                if ~isempty(trail_after_handle) && ishandle(trail_after_handle)
                    delete(trail_after_handle);
                end
                trail_after_handle = plot(X_results, H_results, 'k:', 'LineWidth', 3, 'DisplayName', '潜在升空轨迹');
            end
            legend('Location', 'best');
        end

        % -------- 更新图二：系留绳张力分布 --------
        figure(fig2);

        % 清除上一次查询结果
        if ~isempty(T_query) && ishandle(T_query)
            delete(T_query);
        end
        
        % 获取查询点绳段数据
        L_def_seg_q = tet_data(:,4);                      % 各段形变后长度
        cum_L_from_CP_q = cumsum([0; L_def_seg_q]);       % 距主系点的长度
        L_from_ground_q = L_tet_def - cum_L_from_CP_q;    % 距地面的长度
        T_nodes_q = [tet_data(:,1); F_ground];            % 绳结点张力

        % 查询点绳段张力绘制
        T_query = plot(L_from_ground_q, T_nodes_q, 'b-', 'LineWidth', 3, 'DisplayName', '查询点绳内张力');

        % ====== 更新图二：主系点拉力 ======
        idx_T = find(L_tet_def_results <= L_tet_def, 1, 'last');
        if ~isempty(idx_T)
            if idx_T >= 1
                set(T_CP_handle, 'XData', L_tet_def_results(1:idx_T), 'YData', T_CP_results(1:idx_T));
                if ~isempty(T_CP_after_handle) && ishandle(T_CP_after_handle)
                    delete(T_CP_after_handle);
                end
                if idx_T < length(L_tet_def_results)
                    T_CP_after_handle = plot(L_tet_def_results(idx_T+1:end), T_CP_results(idx_T+1:end), ...
                        'k:', 'LineWidth', 3, 'DisplayName', '主系点拉力');
                else
                    T_CP_after_handle = [];
                end
            else
                set(T_CP_handle, 'XData', [], 'YData', []);
                if ~isempty(T_CP_after_handle) && ishandle(T_CP_after_handle)
                    delete(T_CP_after_handle);
                end
                T_CP_after_handle = plot(L_tet_def_results, T_CP_results, 'k:', 'LineWidth', 3, ...
                    'DisplayName', '潜在主系点拉力');
            end
            legend('Location', 'best');
        end
    end
    fprintf('======== 定高查询结束！========\n\n');
end
fprintf('程序已退出。\n');

%% ======== 辅助函数 ========
% -------- 定高求解函数 --------
function [alpha, T_CP, theta, drift, tet_data, valid, ...
    F_b, W_ast, F_l, F_d, M_CG, L_tet_ori, L_tet_def, W_tet, ...
    g, T, P, rho_air, rho_He, mu, v, z_CP, x_cum, F_ground, ...
    dM_dv, dM_da] = SolveAtHeight(...
    h, L_ast, V_ast, S_ref, m_enve, m_load, x_CB, y_CB, ...
    d_tet, A_tet, rho_tet, L_seg, E, x_CP, y_CP, ...
    v_ref, h_ref, h_ground, alpha_ws, ...
    R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, alpha0)

    % 初始化返回值
    alpha = 0;  F_l = 0;  F_d = 0;  M_CG = 0;  dM_dv = 0;  dM_da = 0;  drift = 0;
    tet_data = [];  T_CP = 0;  theta = 0;  L_tet_ori = 0;  L_tet_def = 0;  W_tet = 0;
    z_CP = 0;  x_cum = 0;  F_ground = 0;  valid = false;
    
    % 计算浮空器某高度的环境参数
    [g, T, P, rho_air, rho_He, mu, v] = ISA_params( ...
        h, R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, v_ref, h_ref, alpha_ws);
    
    % 浮空器浮力与总重力
    m_He = rho_He * V_ast;                   % 氦气质量
    W_ast = (m_enve + m_load + m_He) * g;    % 浮空器总重力
    F_b = rho_air * V_ast * g;               % 浮空器浮力
    if F_b < W_ast                           % 检查浮空器浮力是否小于重力
        return;
    end
    
    % 使用fsolve求解力矩方程
    if nargin < 23 || isempty(alpha0)
        alpha0 = 0;    % 若未提供初始猜测，使用默认值
    end

    options = optimoptions('fsolve', 'Display', 'off');    % fsolve求解设置
    try
        [alpha, fval, exitflag] = fsolve(@(a) Moment_eq(a, ...
            rho_air, v, S_ref, L_ast, F_b, W_ast, x_CB, y_CB, x_CP, y_CP), alpha0, options);
    catch
        return;
    end
    if exitflag <= 0 || abs(fval) > 1e-3    % 检查收敛性和残差值
        return;
    end

    % 调用气动参数计算拉力和夹角
    [F_l, F_d, M_CG, Cl, Cd, Cm, dCl, dCd, dCm, q] = Aero_params(rho_air, v, S_ref, L_ast, alpha);
    T_CP = sqrt(F_d^2 + (F_l + F_b - W_ast)^2);    % 系留绳拉力
    theta = atan2(F_d, F_l + F_b - W_ast);         % 拉力夹角
    if T_CP < 0 || theta < 0 || theta > pi/2       % 检查拉力、夹角是否有物理意义
        return;
    end

    % 主系点坐标转换（体轴系到风轴系）
    x_CP_w = x_CP*cos(alpha) - y_CP*sin(alpha);
    y_CP_w = x_CP*sin(alpha) + y_CP*cos(alpha);

    % 迎风低头判据（动压静稳定性导数）
    dF_l_dv = rho_air * v * S_ref * Cl;
    dF_d_dv = rho_air * v * S_ref * Cd;
    dM_CG_dv = rho_air * v * S_ref * L_ast * Cm;
    dM_dv = dM_CG_dv - x_CP_w * dF_l_dv - y_CP_w * dF_d_dv;    % 动压静稳定性导数

    % 迎角静稳定性判据
    dF_l_da = q * S_ref * dCl;
    dF_d_da = q * S_ref * dCd;
    dM_CG_da = q * S_ref * L_ast * dCm;
    dx_CP_w_da = -x_CP*sin(alpha) - y_CP*cos(alpha);
    dy_CP_w_da =  x_CP*cos(alpha) - y_CP*sin(alpha);
    dx_CB_w_da = -x_CB*sin(alpha) - y_CB*cos(alpha);
    dM_da = dM_CG_da + F_b * dx_CB_w_da ...                    % 迎角静稳定性导数
            - (F_l + F_b - W_ast) * dx_CP_w_da - dF_l_da * x_CP_w ...
            - F_d * dy_CP_w_da - dF_d_da * y_CP_w;

    % 初始化系留绳求解变量
    z_CP = h + y_CP * cos(alpha) + x_CP * sin(alpha);    % 主系点海拔高度
    F_seg = T_CP;                                        % 绳段张力
    F_next = 0;                                          % 下一段张力
    theta_seg = theta;                                   % 张力夹角
    z_up = z_CP;                                         % 绳段上端点高度
    h_cum = 0;                                           % 系留绳竖直投影
    x_cum = 0;                                           % 系留绳水平投影
    L_seg_def_cum = 0;                                   % 绳段累积形变后长度
    W_seg_cum = 0;                                       % 绳段累积重力
    tet_data = zeros(0,7);                               % 系留绳数据存储

    % 自上而下求解系留绳
    while true
        % 绳段形变
        L_seg_def = L_seg + (F_seg * L_seg) / (E * A_tet);    % 绳段形变后长度
        d_seg_def = d_tet * sqrt(L_seg / L_seg_def);          % 绳段形变后直径
        stress_seg = F_seg / (pi/4 * d_seg_def^2) / 1e6;      % 绳段应力 (MPa)
        m_seg = rho_tet * L_seg * pi/4 * d_tet^2;             % 绳段质量
    
        % 判断下端点是否触地
        z_down = z_up - L_seg_def * cos(theta_seg);
        if z_down < h_ground
            drift = x_cum + x_CP * cos(alpha) - y_CP * sin(alpha);   % 浮空器漂移距离
            n_seg = size(tet_data, 1);                               % 已添加的绳段数量
            L_tet_ori = n_seg * L_seg;                               % 系留绳原长
            L_tet_def = L_seg_def_cum;                               % 系留绳形变后长度
            W_tet = W_seg_cum;                                       % 系留绳重力
            F_ground = F_next;                                       % 地面锚点拉力
            valid = true;
            break;
        end
    
        % 未触地则计算绳段中点环境参数
        z_mid = z_up - 0.5 * L_seg_def * cos(theta_seg);
        [g_seg, ~, ~, rho_air_seg, ~, mu_seg, v_seg] = ISA_params( ...
            z_mid, R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, v_ref, h_ref, alpha_ws);
    
        % 绳段法向阻力
        v_n = v_seg * cos(theta_seg);
        D_n = 0;
        if abs(v_n) > 0.01
            Re_n = rho_air_seg * abs(v_n) * d_seg_def / mu_seg;
            if     Re_n < 30;     Cd_n = 1.45 + 8.55 * Re_n^(-0.4);
            elseif Re_n < 100;    Cd_n = 1 + 4 * Re_n^(-0.5);
            elseif Re_n < 1000;   Cd_n = 2.25 - 0.45 * log10(Re_n);
            elseif Re_n < 4000;   Cd_n = 0.9;
            elseif Re_n < 15000;  Cd_n = 1.05 + 0.054 * log10(Re_n);
            elseif Re_n < 150000; Cd_n = 1.21;
            else;                 Cd_n = 0.3;
            end
            D_n = 0.5 * rho_air_seg * v_n^2 * L_seg_def * d_seg_def * Cd_n;
        end
    
        % 绳段切向阻力
        v_t = v_seg * sin(theta_seg);
        D_t = 0;
        if abs(v_t) > 0.001
            Re_t = rho_air_seg * abs(v_t) * d_seg_def / mu_seg;
            if Re_t < 100.55; Cd_t = 1.88 / Re_t^0.74;
            else;             Cd_t = 0.062;
            end
            D_t = 0.5 * rho_air_seg * v_t^2 * L_seg_def * d_seg_def * Cd_t;
        end
    
        % 微元段平衡方程
        F_horiz = F_seg * sin(theta_seg) + D_n * cos(theta_seg) + D_t * sin(theta_seg);
        F_vert = F_seg * cos(theta_seg) - m_seg * g_seg - D_n * sin(theta_seg) + D_t * cos(theta_seg);
        if F_vert < 0    % 垂直力为负表示张力向上，系留绳松弛，浮空器无法升空
            return;
        end
    
        F_next = sqrt(F_horiz^2 + F_vert^2);    % 下一段张力
        theta_next = atan(F_horiz / F_vert);    % 下一段张力与竖直方向夹角
    
        % 记录当前段数据（7列：张力，角度，下端点高度，形变后长度，应力，法向阻力，切向阻力）
        tet_data(end+1, :) = [F_seg, theta_seg, z_down, L_seg_def, stress_seg, D_n, D_t];
    
        % 微元段累积
        h_cum = h_cum + L_seg_def * cos(theta_seg);    % 系留绳竖直投影
        x_cum = x_cum + L_seg_def * sin(theta_seg);    % 系留绳水平投影
        L_seg_def_cum = L_seg_def_cum + L_seg_def;     % 系留绳形变后长度
        W_seg_cum = W_seg_cum + m_seg * g_seg;         % 系留绳总重力
    
        % 更新下一段值
        F_seg = F_next;
        theta_seg = theta_next;
        z_up = z_CP - h_cum;
    end
end

% -------- 对主系点的力矩残差函数 --------
function res = Moment_eq(alpha, rho_air, v, S_ref, L_ast, F_b, W_ast, x_CB, y_CB, x_CP, y_CP)
    % 体轴系到风轴系的坐标转换
    x_CP_w = x_CP*cos(alpha) - y_CP*sin(alpha);
    y_CP_w = x_CP*sin(alpha) + y_CP*cos(alpha);
    x_CB_w = x_CB*cos(alpha) - y_CB*sin(alpha);

    [F_l, F_d, M_CG] = Aero_params(rho_air, v, S_ref, L_ast, alpha);      % 调用气动力计算函数

    res = M_CG + F_b*x_CB_w - (F_l + F_b - W_ast)*x_CP_w - F_d*y_CP_w;    % 对主系点的力矩平衡方程
end

% -------- 浮空器气动力计算函数--------
function [F_l, F_d, M_CG, Cl, Cd, Cm, dCl, dCd, dCm, q] = Aero_params(rho_air, v, S_ref, L_ast, alpha)
    % ★多项式系数★
    Cl_coeff = [1.273, 0];    % 从左往右依次为高次到低次，最右边为常数项
    Cd_coeff = [1.4235, -0.0134, 0.0356];
    Cm_coeff = [11.333, -0.2994, -6.0303, 0.0329, 0.4971, -0.0012];

    Cl = polyval(Cl_coeff, alpha);    % 升力系数
    Cd = polyval(Cd_coeff, alpha);    % 阻力系数
    Cm = polyval(Cm_coeff, alpha);    % 俯仰力矩系数

    % 气动系数求导
    dCl = polyval(polyder(Cl_coeff), alpha);
    dCd = polyval(polyder(Cd_coeff), alpha);
    dCm = polyval(polyder(Cm_coeff), alpha);

    q = 0.5 * rho_air * v^2;          % 动压
    F_l =  q * S_ref * Cl;            % 气动升力
    F_d =  q * S_ref * Cd;            % 气动阻力
    M_CG = q * S_ref * Cm * L_ast;    % 俯仰力矩【使浮空器抬头为正】
end

% -------- 大气环境参数计算函数 --------
function [g, T, P, rho_air, rho_He, mu, v] = ISA_params( ...
    h, R_earth, g0, P0, T0, gamma, R_air, R_He, mu0, S, v_ref, h_ref, alpha_ws)
    % ISA模型
    g = g0 * (R_earth / (R_earth + h))^2;            % 重力加速度 (m/s²) 
    T = T0 - gamma * h;                              % 温度 (K)【对流层温度变化规律】
    P = P0 * (T / T0)^(g / (R_air * gamma));         % 压强 (Pa)【由静力学方程与理想气体状态方程推导】
    rho_air = P / (R_air * T);                       % 空气密度 (kg/m³)【理想气体状态方程；海平面值为1.225】
    rho_He = P / (R_He * T);                         % 氦气密度 (kg/m³)
    mu = mu0 * (T / T0)^1.5 * (T0 + S) / (T + S);    % 动力粘度 (Pa·s) 【萨瑟兰公式】
    v = v_ref * (h / h_ref)^alpha_ws;                % 风速 (m/s) 【幂律风廓线】
end