clc; clear; close all;
# SOLUCAO EM PROGRESSO - CODIGO COM ERRO OU INCERTEZA DA SOLUCAO
% Resolução da Equação de Van der Pol com o Método de Newton

% O problema é descrito pela equação diferencial não linear:
% u'' - μ(u^2 - 1)u' + u = 0, com μ = 1/2, u(0) = 0 e u'(2) = 1.

% Este código resolve o sistema não linear utilizando o método de Newton para encontrar uma solução numérica.
% Este código usa o método de Newton para resolver a equação de Van der Pol de forma numérica.
% A solução é obtida por meio de uma discretização das equações diferenciais de primeira ordem e a aplicação iterativa do método de Newton.
% A convergência é monitorada através de uma condição de parada, e os resultados são visualizados em um gráfico.
% Passos seguidos:

% Definição dos parâmetros e discretização do intervalo [0, 2].
% Construção do sistema não linear utilizando o método de diferenças finitas.
% Implementação do método de Newton para resolver o sistema não linear.
% Exibição do gráfico da solução aproximada para u(t).

% Definição dos parâmetros iniciais
mu = 1/2;         % Parâmetro μ da equação de Van der Pol
t0 = 0; t_end = 2; % Intervalo [0, 2]
N = 50;            % Número de passos (subintervalos)
h = (t_end - t0) / (N + 1); % Passo h, com N divisões
t = linspace(t0, t_end, N+2); % Vetor de pontos de malha de tempo

% Condições iniciais: u(0) = 0 e u'(2) = 1
U = zeros(N, 1);  % Aproximação para u(t) nos pontos internos
V = ones(N, 1);   % Aproximação para u'(t) nos pontos internos

% Método de Newton
tolerance = 1e-4;  % Tolerância para a convergência
max_iter = 100;    % Número máximo de iterações para o método de Newton

% Vetor para armazenar o erro a cada iteração
erro_iteracao = zeros(max_iter, 1);

% Iteração do método de Newton
for iter = 1:max_iter
    F = zeros(2*N, 1); % Vetor das equações no sistema não linear
    J = zeros(2*N, 2*N); % Matriz Jacobiana do sistema

    % Montando o sistema não linear F e sua Jacobiana J
    for i = 1:N-1
        % Equação 1: U_{i+1} - U_i = h * V_i
        F(2*i-1) = U(i+1) - U(i) - h * V(i);

        % Equação 2: V_{i+1} - V_i = h * (mu * (U_i^2 - 1) * V_i - U_i)
        F(2*i) = V(i+1) - V(i) - h * (mu * (U(i)^2 - 1) * V(i) - U(i));

        % Derivadas parciais da Equação 1
        J(2*i-1, 2*i-1) = 1;  % dF1/dU_{i+1}
        J(2*i-1, 2*i) = -h;   % dF1/dV_{i+1}

        % Derivadas parciais da Equação 2
        J(2*i, 2*i-1) = -h * (2 * mu * U(i) * V(i) - 1);  % dF2/dU_{i+1}
        J(2*i, 2*i) = 1 + h * mu * (U(i)^2 - 1);          % dF2/dV_{i+1}
    end

    % Condição de contorno u'(2) = 1
    F(2*N-1) = V(N) - 1;
    J(2*N-1, 2*N-1) = 1;  % Derivada da condição de contorno em relação a V(N)

    % Verificação da singularidade da matriz Jacobiana
    if cond(J) > 1e12
        error('Matriz Jacobiana é singular ou quase singular.');
    end

    % Solução do sistema linear J * delta = -F (método de Newton)
    delta = -J \ F;

    % Atualizando as soluções U e V
    U = U + delta(1:2:end);
    V = V + delta(2:2:end);

    % Calculando o erro entre iterações consecutivas
    if iter > 1
        erro_iteracao(iter) = norm(delta, inf);
    end

    % Critério de parada: verifica se a norma do vetor de alterações é suficientemente pequena
    if norm(delta, inf) < tolerance
        fprintf('Convergência alcançada após %d iterações.\n', iter);
        erro_iteracao = erro_iteracao(1:iter);  % Cortando os valores extras de erro
        break;
    end
end

if iter == max_iter
    fprintf('Número máximo de iterações alcançado.\n');
end

% Plotando os resultados
figure;

subplot(2, 1, 1);
plot(t(2:end-1), U, '-o', 'LineWidth', 2);
xlabel('t');
ylabel('u(t)');
title('Aproximação de u(t) usando o Método de Newton');
grid on;
legend('u(t)', 'Location', 'northeast');

subplot(2, 1, 2);
plot(1:iter, erro_iteracao(1:iter), '-x', 'LineWidth', 2);
xlabel('Iterações');
ylabel('Erro');
title('Erro entre iterações');
grid on;
legend('Erro', 'Location', 'northeast');
