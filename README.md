# Older Kidney Transplant Microsimulation Model

## Background
In the United States, those 65 and older make up a large and growing proportion of End-Stage Kidney Disease (ESKD) patients (~50% of all ESKD patients). The primary treatment options for ESKD are dialysis or kidney transplantation. Kidney transplantation is the preferred method with benefits for life expectancy and quality of life, but donor kidneys are a scare resource, which necessitates the transplant waitlist. Unfortunately, older candidates receive a disproportionately low share of transplants. 

A paradox exists where there are long wait times for transplantation, but each year many usable donor kidneys are discarded. This project posits the question, "What if we utilized the often discarded, sub-optimal kidneys for older transplant candidates who are less likely to receive a transplant otherwise?"

## Methods
We utilized Scientific Registry of Transplant Recipient (SRTR) patient-level data to develop risk prediction equations for key transplant outcomes:
- Waitlist Outcomes
  - Deceased Donor Transplant
  - Living Donor Transplant
  - Waitlist Mortality
  - Other Waitlist Removal
- 30-Day Outcomes
  - Graft Success (No Complications)
  - Delayed Gradt Function (Graft Functions by Day 30)
  - Graft Loss
  - Mortality
- Post-Transplant Outcomes
  - Death with Functioning Graft
  - Graft Loss
  - Death Same Day as Graft Loss
  - Death after Graft Loss

We performed calibration to ensure that our model simultaneously matched key observed target outcomes and to reflect uncertainty from these empirical targets via our model parameters to our model-predicted outcomes. Our targets were derived from the observed outcomes of the data held out from equation development. We varied the equation coefficients probabilistically while preserving the correlation structure. This procedure was repeated 4,000,000 times, resulting in 4,000,000 parameter sets. We then utilize acceptance sampling with the goal of generating good coverage of the width of our targets by retaining the top 100 best fitting parameter sets.
