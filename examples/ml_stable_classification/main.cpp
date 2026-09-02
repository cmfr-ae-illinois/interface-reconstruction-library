#include <vtkCellCenters.h>
#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/data_gen.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;
struct SimulationFolders {
    std::string name;
    std::string dataDirectory;
    std::string plicDirectory;
    int downsample_factor;
};

static std::string read_next_nonempty_line(std::ifstream& input)
{
    std::string line;

    while (std::getline(input, line)) {
        // Remove possible Windows carriage return
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (!line.empty()) {
            return line;
        }
    }

    throw std::runtime_error(
        "Unexpected end of classification input file."
    );
}


static void loadClassificationInput(
    const std::string& input_file,
    std::string& dataset_path,
    std::vector<SimulationFolders>& simulations)
{
    std::ifstream input(input_file);

    if (!input) {
        throw std::runtime_error(
            "Could not open classification input file: "
            + input_file
        );
    }

    // Dataset
    dataset_path = read_next_nonempty_line(input);

    // Number of simulations
    const std::string no_simulations_line =
        read_next_nonempty_line(input);

    int no_simulations;

    try {
        no_simulations = std::stoi(no_simulations_line);
    }
    catch (...) {
        throw std::runtime_error(
            "Invalid number of simulations in input file: "
            + no_simulations_line
        );
    }

    if (no_simulations < 1) {
        throw std::runtime_error(
            "Number of simulations must be at least 1."
        );
    }

    // Simulations

    simulations.clear();
    simulations.reserve(no_simulations);

    for (int i = 0; i < no_simulations; ++i) {
        SimulationFolders simulation;

        simulation.name =
            read_next_nonempty_line(input);

        simulation.dataDirectory =
            read_next_nonempty_line(input);

        simulation.plicDirectory =
            read_next_nonempty_line(input);

        const std::string downsample_line =
            read_next_nonempty_line(input);

        try {
            simulation.downsample_factor =
                std::stoi(downsample_line);
        }
        catch (...) {
            throw std::runtime_error(
                "Invalid downsample factor for simulation "
                + simulation.name
                + ": "
                + downsample_line
            );
        }

        if (simulation.downsample_factor < 1) {
            throw std::runtime_error(
                "Downsample factor for simulation "
                + simulation.name
                + " must be at least 1."
            );
        }

        simulations.push_back(
            std::move(simulation)
        );
    }
}

// Fraction of equal entries between two prediction vectors
static double agreement_fraction(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Prediction vectors have different sizes; align by cellId instead.");
    }
    if (a.empty()) return 0.0;

    size_t same = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        same += (a[i] == b[i]);
    }
    return static_cast<double>(same) / static_cast<double>(a.size());
}

static int pick_most_agreeing_model(
    const std::vector<std::vector<int>>& predictions,
    double* out_mean_agreement = nullptr)
{
    const int R = static_cast<int>(predictions.size());

    if (R == 0) {
        if (out_mean_agreement) {
            *out_mean_agreement = 0.0;
        }
        return -1;
    }

    if (R == 1) {
        if (out_mean_agreement) {
            *out_mean_agreement = 1.0;
        }
        return 1;  // run 1
    }

    const size_t N = predictions[0].size();

    for (int run = 2; run <= R; ++run) {
        if (predictions[run - 1].size() != N) {
            throw std::runtime_error(
                "Prediction vectors differ in length!"
            );
        }
    }

    int most_agreeing_run = 1;
    double best_score = -1.0;

    // Run numbers are 1, ..., R
    for (int run = 1; run <= R; ++run) {
        double sum_agree = 0.0;

        for (int other_run = 1; other_run <= R; ++other_run) {
            if (other_run == run) {
                continue;
            }

            sum_agree += agreement_fraction(
                predictions[run - 1],
                predictions[other_run - 1]
            );
        }

        const double mean_agree =
            sum_agree / static_cast<double>(R - 1);

        if (mean_agree > best_score) {
            best_score = mean_agree;
            most_agreeing_run = run;
        }
    }

    if (out_mean_agreement) {
        *out_mean_agreement = best_score;
    }

    return most_agreeing_run;
}

// Compute and print mean per-cell instability
static double compute_mean_instability(const std::vector<std::vector<int>>& predictions) {
    const int R = static_cast<int>(predictions.size());
    if (R == 0) return 0.0;

    const size_t N = predictions[0].size();
    for (int r = 1; r < R; ++r) {
        if (predictions[r].size() != N) {
            throw std::runtime_error("Prediction vectors differ in length!");
        }
    }

    double sum_u = 0.0;
    for (size_t t = 0; t < N; ++t) {
        std::array<int,6> counts{0,0,0,0,0,0};
        for (int r = 0; r < R; ++r) {
            int c = predictions[r][t];
            if (c < 0 || c >= 6) continue; // should not happen
            counts[c]++;
        }
        int max_count = *std::max_element(counts.begin(), counts.end());
        double u = 1.0 - static_cast<double>(max_count) / static_cast<double>(R);
        sum_u += u;
    }

    return (N > 0) ? (sum_u / static_cast<double>(N)) : 0.0;
}

int main (int argc, char* argv[]) {

    const int no_runs = 10; // Number of independent runs for stable classification, number of models to train and evaluate
    double pdistribution_step = 0.0; // Not needed for stable classification, but required by classify_simulation

    const std::string default_input_file =
    "../../../examples/ml_stable_classification/stable_classification_input.txt";

    std::string input_file;

    if (argc >= 2) {
        input_file = argv[1];

        std::cout
            << "Using classification input file: "
            << input_file
            << "\n";
    }
    else {
        input_file = default_input_file;

        std::cout
            << "No classification input file provided. "
            << "Using default:\n"
            << input_file
            << "\n";
    }

    std::string dataset_path;
    std::vector<SimulationFolders> simulations;

    loadClassificationInput(
        input_file,
        dataset_path,
        simulations
    );

    std::cout
        << "\n========================================\n"
        << "Classification input\n"
        << "========================================\n";

    std::cout
        << "Dataset: "
        << dataset_path
        << "\n";

    std::cout
        << "Number of simulations: "
        << simulations.size()
        << "\n";

    for (std::size_t i = 0; i < simulations.size(); ++i) {
        const auto& simulation = simulations[i];

        std::cout
            << "\nSimulation " << (i + 1) << ":\n"
            << "  Name       : " << simulation.name << "\n"
            << "  Data       : " << simulation.dataDirectory << "\n"
            << "  PLIC       : " << simulation.plicDirectory << "\n"
            << "  Downsample : " << simulation.downsample_factor << "\n";
    }

    std::cout
        << "========================================\n";

    const auto now = std::chrono::system_clock::now();
    const std::time_t currentTime =
        std::chrono::system_clock::to_time_t(now);

    std::tm localTime{};

#if defined(_WIN32)
    localtime_s(&localTime, &currentTime);
#else
    localtime_r(&currentTime, &localTime);
#endif

    char timestamp[64];

    std::snprintf(
        timestamp,
        sizeof(timestamp),
        "%04d-%02d-%02d_%02d%02d%02d",
        localTime.tm_year + 1900,
        localTime.tm_mon + 1,
        localTime.tm_mday,
        localTime.tm_hour,
        localTime.tm_min,
        localTime.tm_sec
    );

    // Create a distinct output directory
    fs::path experiment_dir =
        fs::path("stable_run_models") / timestamp;

    fs::create_directories(experiment_dir);

    // Stored results

    // predictions[r] contains every classified cell from every timestep of every simulation, in deterministic order
    std::vector<std::vector<int>> predictions;
    predictions.reserve(no_runs);

    // Directory containing each trained model
    std::vector<fs::path> run_dirs;
    run_dirs.reserve(no_runs);

    // Total classified-cell count for each simulation in the first run
    std::vector<std::size_t> reference_simulation_cell_counts;

    // Per-timestep classified-cell counts for every simulation
    // First index: simulation
    // Second index: timestep
    std::vector<std::vector<std::size_t>>
        reference_timestep_cell_counts;

    // Train and evaluate all model runs
    for (int run = 1; run < no_runs+1; ++run) {
        std::cout
            << "\n==================================================\n"
            << "Stable classification run "
            << (run)
            << " / "
            << no_runs
            << "\n"
            << "==================================================\n";

        fs::path run_dir = experiment_dir / ("run_" + std::to_string(run));

        fs::create_directories(run_dir);
        run_dirs.push_back(run_dir);

        // Create a fresh classifier for this run.
        IRL::MLClassifier ml;

        ml.loadDataset(dataset_path);

        ml.preprocess_data();

        // A deterministic run-specific seed could be set here:
        //
        ml.setSeed(1234 + run);

        ml.trainModel();

        // Contains every prediction from every simulation and timestep.
        std::vector<int> savedClassesCombined;

        // Total classified cells for each simulation in this run.
        std::vector<std::size_t> current_simulation_cell_counts;
        current_simulation_cell_counts.reserve(simulations.size());

        // Per-timestep counts for each simulation in this run.
        std::vector<std::vector<std::size_t>>
            current_timestep_cell_counts;

        current_timestep_cell_counts.reserve(simulations.size());

        // -------------------------------------------------------------
        // Classify every timestep of every simulation
        // -------------------------------------------------------------

        for (const auto& simulation : simulations) {
            std::cout
                << "\nClassifying complete simulation: "
                << simulation.name
                << "\n";

            std::vector<int> savedClassesThisSimulation;

            std::vector<std::size_t>
                timestepCellCountsThisSimulation;

            // No classified VTK files are written during the ten model-
            // selection runs. Only predictions are collected.
            IRL::classify_simulation(
                ml,
                simulation.dataDirectory,
                simulation.plicDirectory,
                "", // unused because write_output is false
                &savedClassesThisSimulation,
                simulation.downsample_factor,
                pdistribution_step,
                false, // write_output
                &timestepCellCountsThisSimulation
            );

            const std::size_t simulationCellCount =
                savedClassesThisSimulation.size();

            current_simulation_cell_counts.push_back(
                simulationCellCount
            );

            current_timestep_cell_counts.push_back(
                timestepCellCountsThisSimulation
            );

            std::cout
                << simulation.name
                << " timesteps classified: "
                << timestepCellCountsThisSimulation.size()
                << "\n";

            std::cout
                << simulation.name
                << " total classified interface cells: "
                << simulationCellCount
                << "\n";

            for (std::size_t timestepIndex = 0;
                 timestepIndex <
                     timestepCellCountsThisSimulation.size();
                 ++timestepIndex) {

                std::cout
                    << "  timestep "
                    << IRL::makeTimestepTag(
                           static_cast<int>(timestepIndex + 1))
                    << ": "
                    << timestepCellCountsThisSimulation[timestepIndex]
                    << " classified cells\n";
            }

            // Append this complete simulation to the combined vector.
            savedClassesCombined.insert(
                savedClassesCombined.end(),
                savedClassesThisSimulation.begin(),
                savedClassesThisSimulation.end()
            );
        }

        // -------------------------------------------------------------
        // Verify that all model runs classified identical cells
        // -------------------------------------------------------------

        if (run == 1) {
            reference_simulation_cell_counts =
                current_simulation_cell_counts;

            reference_timestep_cell_counts =
                current_timestep_cell_counts;
        } else {
            if (current_simulation_cell_counts.size() !=
                reference_simulation_cell_counts.size()) {

                throw std::runtime_error(
                    "The number of classified simulations differs "
                    "between model runs."
                );
            }

            for (std::size_t simulationIndex = 0;
                 simulationIndex < simulations.size();
                 ++simulationIndex) {

                if (current_simulation_cell_counts[simulationIndex] !=
                    reference_simulation_cell_counts[simulationIndex]) {

                    throw std::runtime_error(
                        "Total classified-cell count differs between "
                        "runs for simulation "
                        + simulations[simulationIndex].name
                    );
                }

                const auto& currentTimestepCounts =
                    current_timestep_cell_counts[simulationIndex];

                const auto& referenceTimestepCounts =
                    reference_timestep_cell_counts[simulationIndex];

                if (currentTimestepCounts.size() !=
                    referenceTimestepCounts.size()) {

                    throw std::runtime_error(
                        "Number of classified timesteps differs between "
                        "runs for simulation "
                        + simulations[simulationIndex].name
                    );
                }

                for (std::size_t timestepIndex = 0;
                     timestepIndex < referenceTimestepCounts.size();
                     ++timestepIndex) {

                    if (currentTimestepCounts[timestepIndex] !=
                        referenceTimestepCounts[timestepIndex]) {

                        throw std::runtime_error(
                            "Classified-cell count differs between runs "
                            "for simulation "
                            + simulations[simulationIndex].name
                            + ", timestep "
                            + IRL::makeTimestepTag(
                                static_cast<int>(
                                    timestepIndex + 1))
                        );
                    }
                }
            }
        }

        if (run > 1 &&
            savedClassesCombined.size() != predictions[0].size()) {

            throw std::runtime_error(
                "Combined prediction-vector size differs between runs."
            );
        }

        predictions.push_back(
            std::move(savedClassesCombined)
        );

        std::cout
            << "\nCombined classified interface cells for run "
            << run
            << ": "
            << predictions.back().size()
            << "\n";

        // Save this trained model.
        ml.saveModel(
            run_dir.string() + "/",
            false
        );

        // ml goes out of scope here, releasing the loaded dataset.
    }

    // ---------------------------------------------------------------------
    // Measure ensemble instability
    // ---------------------------------------------------------------------

    const double mean_instability =
        compute_mean_instability(predictions);

    std::cout
        << "\nMean per-cell instability: "
        << mean_instability
        << "\n";

    // ---------------------------------------------------------------------
    // Select the most agreeing model

    double most_agreeing_mean_agreement = 0.0;

    const int most_agreeing_run =
        pick_most_agreeing_model(
            predictions,
            &most_agreeing_mean_agreement
        );

    if (most_agreeing_run < 0) {
        throw std::runtime_error(
            "No model runs were available for model selection."
        );
    }

    const fs::path most_agreeing_model_dir =
        run_dirs.at(
            static_cast<std::size_t>(most_agreeing_run-1));

    const fs::path most_agreeing_model_path =
        most_agreeing_model_dir / "ml_model.pt";

    // ---------------------------------------------------------------------
    // Write model-selection summary
    // ---------------------------------------------------------------------

    const fs::path selection_file =
        experiment_dir / "model_selection.txt";

    {
        std::ofstream output(selection_file);

        if (!output) {
            throw std::runtime_error(
                "Failed to open model selection file: "
                + selection_file.string()
            );
        }

        output
            << "dataset_path "
            << dataset_path
            << "\n";

        output
            << "no_runs "
            << no_runs
            << "\n";

        output
            << "simulations_used_for_selection "
            << simulations.size()
            << "\n";

        for (std::size_t simulationIndex = 0;
             simulationIndex < simulations.size();
             ++simulationIndex) {

            const auto& simulation =
                simulations[simulationIndex];

            output
                << "simulation_"
                << simulationIndex
                << "_name "
                << simulation.name
                << "\n";

            output
                << "simulation_"
                << simulationIndex
                << "_data_directory "
                << simulation.dataDirectory
                << "\n";

            output
                << "simulation_"
                << simulationIndex
                << "_plic_directory "
                << simulation.plicDirectory
                << "\n";

            output
                << "simulation_"
                << simulationIndex
                << "_downsample_factor "
                << simulation.downsample_factor
                << "\n";

            output
                << "simulation_"
                << simulationIndex
                << "_number_of_timesteps "
                << reference_timestep_cell_counts[
                       simulationIndex].size()
                << "\n";

            output
                << "simulation_"
                << simulationIndex
                << "_total_classified_cells "
                << reference_simulation_cell_counts[
                       simulationIndex]
                << "\n";

            for (std::size_t timestepIndex = 0;
                 timestepIndex <
                     reference_timestep_cell_counts[
                         simulationIndex].size();
                 ++timestepIndex) {

                output
                    << "simulation_"
                    << simulationIndex
                    << "_timestep_"
                    << IRL::makeTimestepTag(
                           static_cast<int>(
                               timestepIndex + 1))
                    << "_classified_cells "
                    << reference_timestep_cell_counts[
                           simulationIndex][timestepIndex]
                    << "\n";
            }
        }

        output
            << "most_agreeing_model_run "
            << most_agreeing_run
            << "\n";

        output
            << "most_agreeing_model_path "
            << most_agreeing_model_path.string()
            << "\n";

        output
            << "most_agreeing_mean_agreement "
            << most_agreeing_mean_agreement
            << "\n";

        output
            << "mean_per_cell_instability "
            << mean_instability
            << "\n";
    }

    // ---------------------------------------------------------------------
    // Print model-selection result
    // ---------------------------------------------------------------------

    std::cout
        << "\n==================================================\n"
        << "Model selection result\n"
        << "==================================================\n";

    std::cout
        << "Selection file: "
        << selection_file.string()
        << "\n";

    std::cout
        << "Most agreeing run: "
        << most_agreeing_run
        << "\n";

    std::cout
        << "Most agreeing model: "
        << most_agreeing_model_path.string()
        << "\n";

    std::cout
        << "Most agreeing mean agreement: "
        << most_agreeing_mean_agreement
        << "\n";

    std::cout
        << "Mean per-cell instability: "
        << mean_instability
        << "\n";

    for (std::size_t simulationIndex = 0;
         simulationIndex < simulations.size();
         ++simulationIndex) {

        std::cout
            << simulations[simulationIndex].name
            << ": "
            << reference_timestep_cell_counts[
                   simulationIndex].size()
            << " timesteps, "
            << reference_simulation_cell_counts[
                   simulationIndex]
            << " classified cells\n";
    }

    // ---------------------------------------------------------------------
    // Reload selected model and classify complete simulations with output
    // ---------------------------------------------------------------------

    {
        IRL::MLClassifier most_agreeing;

        most_agreeing.loadModel(
            most_agreeing_model_path.string(),
            false
        );

        const fs::path selectedOutputsRoot =
            experiment_dir / "selected_model_outputs";

        for (const auto& simulation : simulations) {
            const fs::path simulationOutputDirectory =
                selectedOutputsRoot / simulation.name;

            std::cout
                << "\nReclassifying complete simulation "
                << simulation.name
                << " with the selected model.\n";

            std::cout
                << "Output directory: "
                << simulationOutputDirectory.string()
                << "\n";

            IRL::classify_simulation(
                most_agreeing,
                simulation.dataDirectory,
                simulation.plicDirectory,
                simulationOutputDirectory.string(),
                nullptr,
                simulation.downsample_factor,
                pdistribution_step,
                true,    // write_output
                nullptr  // no timestep counts needed
            );
        }
    }
}