// Simulate form submission and show results
document.getElementById('drug-form').addEventListener('submit', function (e) {
    e.preventDefault(); // Prevent actual form submission

    // Show loading spinner
    document.getElementById('button-text').classList.add('hidden');
    document.getElementById('loading-spinner').classList.remove('hidden');

    // Simulate a delay for analysis
    setTimeout(() => {
        // Hide loading spinner
        document.getElementById('button-text').classList.remove('hidden');
        document.getElementById('loading-spinner').classList.add('hidden');

        // Show results
        document.getElementById('result-section').classList.remove('hidden');
        document.getElementById('result-text').textContent = "No interactions found between the provided drugs.";
    }, 2000); // 2-second delay for simulation
});

// Optional: Add dark mode toggle
document.getElementById('theme-toggle').addEventListener('click', function () {
    document.body.classList.toggle('dark-mode');
});