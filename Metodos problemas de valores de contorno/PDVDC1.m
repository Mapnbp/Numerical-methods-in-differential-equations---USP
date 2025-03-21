clc; clear; close all;

# Este código usa diferenças finitas para transformar a EDO em um sistema linear.
# e depois resolve o sistema usando A \ b.

# Segue basicamente os sequintes passos:
% Discretiza a equação diferencial 𝑢′′ = 4(𝑢−𝑥) usando diferenças finitas.
% Resolve o sistema linear resultante.
% Compara os resultados com a solução exata u(x).
% Gera um gráfico mostrando o erro para diferentes valores de h

% Definição da função exata
u_exact = @(x) exp(2)*(exp(4)-1)^(-1) * (exp(2*x) - exp(-2*x)) + x;

% Valores de h para comparação
h_values = [1/2, 1/4];

figure; hold on;
title("Aproximação da Solução e Erro")
xlabel("x")
ylabel("u(x)")
colors = ['r', 'b']; % Cores para h=1/2 e h=1/4

for idx = 1:length(h_values)
    h = h_values(idx);
    N = 1/h - 1; % Número de nós internos
    x = linspace(0, 1, N+2)'; % Nós incluindo os extremos

    % Montagem da matriz do sistema linear Au = b
    A = (1/h^2) * (diag(2*ones(N,1)) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1));
    b = 4*(x(2:end-1) - u_exact(x(2:end-1)));

    % Condições de contorno u(0) = 0 e u(1) = 2
    b(1) = b(1) + 0/h^2; % u(0) = 0 já está correto
    b(end) = b(end) + 2/h^2;

    % Resolução do sistema linear
    u_num = A \ b;

    % Adiciona as condições de contorno na solução
    u_num = [0; u_num; 2];

    % Solução exata nos pontos
    u_real = u_exact(x);

    % Erro absoluto
    erro = abs(u_real - u_num);

    % Gráficos
    subplot(1,2,1);
    plot(x, u_num, ['o-' colors(idx)], 'DisplayName', sprintf("h = %0.2f", h));
    hold on;

    subplot(1,2,2);
    plot(x, erro, ['.-' colors(idx)], 'DisplayName', sprintf("Erro para h = %0.2f", h));
    hold on;
end

subplot(1,2,1);
plot(x, u_exact(x), 'k--', 'DisplayName', 'Solução Exata');
legend();
title("Aproximação da Solução");
xlabel("x");
ylabel("u(x)");

subplot(1,2,2);
legend();
title("Erro Absoluto");
xlabel("x");
ylabel("|u_{exato} - u_{aproximado}|");

hold off;

