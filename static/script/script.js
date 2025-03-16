const body = document.querySelector("body")
const toggle = document.querySelector("#toggle")
const sun = document.querySelector(".toggle .bxs-sun")
const moon = document.querySelector(".toggle .bx-moon")


toggle.addEventListener("change", () => {
    body.classList.toggle("dark");
    console.log(sun.src.split('/').slice(-2).join('/'))
    sun.className = sun.className == "bx bxs-sun" ? "bx bx-sun" : "bx bxs-sun"
    moon.className = moon.className == "bx bx-moon" ? "bx bxs-moon" : "bx bx-moon"
    moon.src = moon.src.split('/').slice(-2).join('/') == "resources/Moon_fill_light.svg" ? "resources/Moon_fill.svg" : "resources/Moon_fill_light.svg"
    sun.src = sun.src.split('/').slice(-2).join('/') == "resources/Sun_fill.svg" ? "resources/Sun_fill_light.svg" : "resources/Sun_fill.svg"
})
