clc; clear; close all;

# Este código segue os seguintes passos:

% Discretiza a equação diferencial 𝑢′′ = -(𝑢)^2 - u + ln(x)
% Monta o sistema NÃO linear F(U) = 0 com F e a jacobiana J(U)
% Resolve com newton: U^(k+1) = U^(k) −J^−1 F(U^(k)).
% Testa valores menores de h para analisar a convergência
% Plota a solução para h= 1/2, 1/4, 1/8.

% Configurações
h_values = [1/4, 1/8, 1/16]; % Testar diferentes valores de h
tolerance = 1e-6;  % Tolerância para o método de Newton
max_iter = 10;     % Número máximo de iterações

figure; hold on;
title("Aproximação da Solução para Diferentes Valores de h")
xlabel("x")
ylabel("u(x)")
grid on;

for idx = 1:length(h_values)
    h = h_values(idx);
    N = round((2 - 1) / h) - 1; % Número de pontos internos
    x = linspace(1, 2, N+2)'; % Incluindo os extremos

    % Chute inicial (Interpolação linear entre u(1)=0 e u(2)=ln(2))
    U = (log(2) / (2 - 1)) * (x(2:end-1) - 1);

    for iter = 1:max_iter
        % Vetor F(U)
        F = zeros(N, 1);
        J = zeros(N, N); % Jacobiana

        for i = 1:N
            xi = x(i+1);
            if i == 1
                u_prev = 0; % Condição de contorno u(1) = 0
            else
                u_prev = U(i-1);
            end

            if i == N
                u_next = log(2); % Condição de contorno u(2) = ln(2)
            else
                u_next = U(i+1);
            end

            % Aproximações de diferenças finitas
            u_prime = (u_next - u_prev) / (2 * h); % Derivada central
            u_double_prime = (u_next - 2*U(i) + u_prev) / (h^2); % Segunda derivada

            % Função F(U) do sistema não linear
            F(i) = u_double_prime + u_prime^2 + U(i) - log(xi);

            % Jacobiana J(U)
            if i > 1
                J(i, i-1) = 1/h^2 - (u_prime / h);
            end
            J(i, i) = -2/h^2 + 1;
            if i < N
                J(i, i+1) = 1/h^2 + (u_prime / h);
            end
        end

        % Atualização de Newton
        delta_U = J \ (-F);
        U = U + delta_U;

        % Critério de parada
        if norm(delta_U, inf) < tolerance
            break;
        end
    end

    % Adiciona valores de contorno
    U = [0; U; log(2)];

    % Plotagem da solução aproximada
    plot(x, U, '-o', 'DisplayName', sprintf("h = %.3f", h));
end

legend();
hold off;
