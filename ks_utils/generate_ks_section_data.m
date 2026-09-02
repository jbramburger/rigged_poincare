function cache = generate_ks_section_data(outputFile)
%GENERATE_KS_SECTION_DATA Generate the event-located KSE section archive.
%
% This function performs one uninterrupted ode45 integration of the K=32
% Fourier--Galerkin model at nu=0.0298. It does not retain the full time
% trajectory. Instead, a streaming output function examines consecutive
% accepted solver steps, locates upward crossings of a2=0 by cubic Hermite
% interpolation, and stores only the crossing times and states.
%
% The original fixed-grid calculation produced 23089 section rows after
% discarding its first crossing. To make the event-located archive directly
% comparable, this function stops after 23090 upward crossings, discards the
% first, and verifies that x(2:2:end,:) contains N=11544 states for F=P^2.
%
% OUTPUT
%   ks_psec_0298_x2=0_event.mat, unless another filename is supplied.
%   The file contains x, m, and nu for compatibility with the original
%   analysis, together with a cache structure containing the full crossing
%   states and numerical provenance.

    if nargin < 1 || isempty(outputFile)
        outputFile = 'ks_psec_0298_x2=0_event.mat';
    end

    nu = 0.0298;
    m = 32;

    initialCondition = [ ...
        0;
        0.49744707;
        2.115777;
       -4.6417336;
       -3.915899;
       -5.449332;
        0.9944877;
        0.06847247;
       -0.20083837;
       -0.872909;
       -0.29952765;
       -0.13332021;
        0.06590602;
       -0.02766516;
        zeros(m-14,1)];

    requestedInterval = [0,12000];
    relativeTolerance = 1e-12;
    absoluteTolerance = 1e-12*ones(1,m);
    maximumStep = 0.1;
    interpolationTolerance = 1e-12;
    validationEndTime = 20;
    targetRawCrossings = 23090;
    discardedCrossings = 1;
    reportEvery = 1000;

    rawCrossingTimes = nan(targetRawCrossings,1);
    rawCrossingStates = nan(targetRawCrossings,m);
    rawCrossingSpeeds = nan(targetRawCrossings,1);
    crossingCount = 0;
    previousTime = [];
    previousState = [];
    lastAcceptedTime = requestedInterval(1);

    rootOptions = optimset( ...
        'TolX',interpolationTolerance, ...
        'Display','off');

    solverOptions = odeset( ...
        'RelTol',relativeTolerance, ...
        'AbsTol',absoluteTolerance, ...
        'MaxStep',maximumStep, ...
        'Refine',1, ...
        'OutputFcn',@collect_crossings);

    fprintf('Generating event-located KSE section data.\n')
    fprintf('The full trajectory is not stored in memory.\n')
    fprintf('Target upward crossings before transient removal: %d\n', ...
        targetRawCrossings)

    ode45(@(t,state) ks_rhs(state,nu,m), ...
        requestedInterval,initialCondition,solverOptions);

    assert(crossingCount == targetRawCrossings, ...
        ['The requested integration interval ended before the target ' ...
        'number of upward crossings was reached.'])

    rawCrossingTimes = rawCrossingTimes(1:crossingCount);
    rawCrossingStates = rawCrossingStates(1:crossingCount,:);
    rawCrossingSpeeds = rawCrossingSpeeds(1:crossingCount);

    retainedIndices = (discardedCrossings+1):crossingCount;
    retainedCrossingTimes = rawCrossingTimes(retainedIndices);
    retainedCrossingStates = rawCrossingStates(retainedIndices,:);
    retainedCrossingSpeeds = rawCrossingSpeeds(retainedIndices);

    sectionCoordinates = [1,3:m];
    x = retainedCrossingStates(:,sectionCoordinates);
    secondReturnData = x(2:2:end,:);

    assert(size(x,1) == 23089, ...
        'The event-located first-return archive must contain 23089 rows.')
    assert(size(secondReturnData,1) == 11544, ...
        'The retained F=P^2 data must contain N=11544 rows.')
    assert(all(retainedCrossingSpeeds > 0), ...
        'A retained crossing violates da2/dt>0.')

    [alternationRate,oddPurity,evenPurity] = ...
        component_alternation(x);
    assert(alternationRate > 0.99 && ...
        oddPurity > 0.99 && evenPurity > 0.99, ...
        'The first-return sequence does not alternate between components.')

    [validationCount,validationTimeError,validationStateError] = ...
        validate_interpolation( ...
        initialCondition,nu,m,validationEndTime,solverOptions, ...
        rawCrossingTimes,rawCrossingStates);

    cache = struct();
    cache.schemaVersion = 1;
    cache.createdOn = datestr(now,30);
    cache.system = 'Kuramoto--Sivashinsky Fourier--Galerkin model';
    cache.nu = nu;
    cache.modes = m;
    cache.initialCondition = initialCondition;
    cache.requestedIntegrationInterval = requestedInterval;
    cache.actualIntegrationInterval = ...
        [requestedInterval(1),rawCrossingTimes(end)];
    cache.lastAcceptedSolverTime = lastAcceptedTime;
    cache.solver = 'ode45';
    cache.relativeTolerance = relativeTolerance;
    cache.absoluteTolerance = absoluteTolerance;
    cache.maximumStep = maximumStep;
    cache.sectionCoordinate = 'a2';
    cache.sectionValue = 0;
    cache.eventDirection = +1;
    cache.crossingRule = 'a2=0 with da2/dt>0';
    cache.eventLocationRule = [ ...
        'Upward sign changes are detected between consecutive accepted ' ...
        'ode45 steps. The crossing is located by cubic Hermite ' ...
        'interpolation using the endpoint states and vector field.'];
    cache.interpolationTolerance = interpolationTolerance;
    cache.validationInterval = [0,validationEndTime];
    cache.validationCrossingCount = validationCount;
    cache.maximumValidationTimeError = validationTimeError;
    cache.maximumValidationStateError = validationStateError;
    cache.targetRawCrossings = targetRawCrossings;
    cache.discardedInitialCrossings = discardedCrossings;
    cache.transientRule = 'Discard the first upward crossing.';
    cache.rawCrossingTimes = rawCrossingTimes;
    cache.rawCrossingStates = rawCrossingStates;
    cache.rawCrossingSpeeds = rawCrossingSpeeds;
    cache.retainedCrossingTimes = retainedCrossingTimes;
    cache.retainedCrossingStates = retainedCrossingStates;
    cache.retainedCrossingSpeeds = retainedCrossingSpeeds;
    cache.sectionData = x;
    cache.firstReturnCount = size(x,1);
    cache.secondReturnCount = size(secondReturnData,1);
    cache.retainedParityForSecondReturn = 2;
    cache.componentAlternationRate = alternationRate;
    cache.oddComponentPurity = oddPurity;
    cache.evenComponentPurity = evenPurity;
    cache.maximumSectionResidual = ...
        max(abs(retainedCrossingStates(:,2)));
    cache.minimumCrossingSpeed = min(retainedCrossingSpeeds);

    save(outputFile,'x','m','nu','cache','-v7')

    fprintf('Saved %s\n',outputFile)
    fprintf('First-return section rows: %d\n',size(x,1))
    fprintf('Second-return rows N: %d\n',size(secondReturnData,1))
    fprintf('Actual integration interval: [%.16g, %.16g]\n', ...
        cache.actualIntegrationInterval)
    fprintf('Maximum |a2| at a retained crossing: %.6e\n', ...
        cache.maximumSectionResidual)
    fprintf('Minimum retained da2/dt: %.16g\n', ...
        cache.minimumCrossingSpeed)
    fprintf('Empirical component alternation rate: %.10f\n', ...
        alternationRate)
    fprintf(['Interpolation check against ode45 Events over [0,%.1f]: ' ...
        '%d crossings\n'],validationEndTime,validationCount)
    fprintf('Maximum validation time error: %.6e\n',validationTimeError)
    fprintf('Maximum validation state error: %.6e\n',validationStateError)

    function stop = collect_crossings(times,states,flag)
        stop = 0;

        if strcmp(flag,'init')
            previousTime = times(1);
            previousState = states(:,1);
            lastAcceptedTime = previousTime;
            return
        elseif strcmp(flag,'done')
            return
        end

        for outputIndex = 1:numel(times)
            currentTime = times(outputIndex);
            currentState = states(:,outputIndex);

            if currentTime <= previousTime
                continue
            end

            if previousState(2) < 0 && currentState(2) >= 0
                stepLength = currentTime-previousTime;
                previousDerivative = ks_rhs(previousState,nu,m);
                currentDerivative = ks_rhs(currentState,nu,m);

                rootFraction = fzero( ...
                    @(fraction) hermite_a2( ...
                    fraction,previousState,currentState, ...
                    previousDerivative,currentDerivative,stepLength), ...
                    [0,1],rootOptions);

                crossingState = hermite_state( ...
                    rootFraction,previousState,currentState, ...
                    previousDerivative,currentDerivative,stepLength);
                crossingTime = previousTime+rootFraction*stepLength;
                crossingDerivative = ks_rhs(crossingState,nu,m);

                if crossingDerivative(2) > 0
                    crossingCount = crossingCount+1;
                    rawCrossingTimes(crossingCount) = crossingTime;
                    rawCrossingStates(crossingCount,:) = crossingState.';
                    rawCrossingSpeeds(crossingCount) = crossingDerivative(2);

                    if mod(crossingCount,reportEvery) == 0
                        fprintf('  crossings = %d, t = %.6f\n', ...
                            crossingCount,crossingTime)
                    end

                    if crossingCount == targetRawCrossings
                        lastAcceptedTime = currentTime;
                        stop = 1;
                        return
                    end
                end
            end

            previousTime = currentTime;
            previousState = currentState;
            lastAcceptedTime = currentTime;
        end
    end
end

function value = hermite_a2( ...
        fraction,leftState,rightState,leftDerivative,rightDerivative,step)
    state = hermite_state( ...
        fraction,leftState,rightState,leftDerivative,rightDerivative,step);
    value = state(2);
end

function state = hermite_state( ...
        fraction,leftState,rightState,leftDerivative,rightDerivative,step)
    fractionSquared = fraction*fraction;
    fractionCubed = fractionSquared*fraction;

    leftValueWeight = 2*fractionCubed-3*fractionSquared+1;
    leftSlopeWeight = fractionCubed-2*fractionSquared+fraction;
    rightValueWeight = -2*fractionCubed+3*fractionSquared;
    rightSlopeWeight = fractionCubed-fractionSquared;

    state = leftValueWeight*leftState ...
        + step*leftSlopeWeight*leftDerivative ...
        + rightValueWeight*rightState ...
        + step*rightSlopeWeight*rightDerivative;
end

function [alternationRate,oddPurity,evenPurity] = ...
        component_alternation(sectionData)
    oddData = sectionData(1:2:end,:);
    evenData = sectionData(2:2:end,:);

    oddCentroid = mean(oddData,1);
    evenCentroid = mean(evenData,1);
    midpoint = 0.5*(oddCentroid+evenCentroid);
    separatingDirection = oddCentroid-evenCentroid;
    scores = (sectionData-midpoint)*separatingDirection.';

    oddPurity = mean(scores(1:2:end) > 0);
    evenPurity = mean(scores(2:2:end) < 0);
    alternationRate = mean(scores(1:end-1).*scores(2:end) < 0);
end

function [crossingCount,timeError,stateError] = validate_interpolation( ...
        initialCondition,nu,m,endTime,streamingOptions, ...
        streamingTimes,streamingStates)
    validationOptions = odeset(streamingOptions, ...
        'OutputFcn',[], ...
        'Events',@positive_a2_event);

    [~,~,eventTimes,eventStates] = ode45( ...
        @(t,state) ks_rhs(state,nu,m), ...
        [0,endTime],initialCondition,validationOptions);

    crossingCount = numel(eventTimes);
    assert(crossingCount > 0, ...
        'The interpolation-validation interval contains no crossings.')
    assert(crossingCount <= numel(streamingTimes), ...
        'The validation run contains more events than the streaming run.')

    timeError = max(abs(eventTimes-streamingTimes(1:crossingCount)));
    stateDifferences = eventStates-streamingStates(1:crossingCount,:);
    stateError = max(vecnorm(stateDifferences,2,2));
end

function [value,isTerminal,direction] = positive_a2_event(~,state)
    value = state(2);
    isTerminal = 0;
    direction = +1;
end

function derivative = ks_rhs(state,nu,modes)
    derivative = zeros(modes,1);

    for k = 1:modes
        derivative(k) = k^2*(1-nu*k^2)*state(k);

        for n = 1:(modes-k)
            derivative(k) = derivative(k) ...
                + 0.5*k*state(n)*state(n+k);
        end

        for j = 1:(k-1)
            derivative(k) = derivative(k) ...
                - 0.25*k*state(j)*state(k-j);
        end
    end
end
