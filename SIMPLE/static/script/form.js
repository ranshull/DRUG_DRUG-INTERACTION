function checkInteractions() {
    const drugInput = document.getElementById("drugInput").value;
    if (!drugInput) {
        alert("Please enter drug names!");
        return;
    }

    fetch("/check_interaction", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ drugs: drugInput })
    })
    .then(response => response.json())
    .then(data => {
        const resultsDiv = document.getElementById("results");
        resultsDiv.innerHTML = "<h3>Results:</h3>";
        data.interactions.forEach(interaction => {
            resultsDiv.innerHTML += `<p>${interaction}</p>`;
        });
    })
    .catch(error => console.error("Error:", error));
}
