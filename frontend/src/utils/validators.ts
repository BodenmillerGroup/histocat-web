export const required = (value) => !!value || "Required";
export const positiveNumber = (value) => (typeof value === "number" && value > 0) || "Should be positive number";
export const nonNegativeNumber = (value) =>
  (typeof value === "number" && value >= 0) || "Should be non-negative number";
export const percentNumber = (value) =>
  (typeof value === "number" && value >= 0 && value <= 100) || "Value should be between 0 and 100";
export const percentFloatNumber = (value) =>
  (typeof value === "number" && value >= 0 && value <= 1) || "Value should be between 0 and 1";

export const email = (value) => {
  const pattern = /^(([^<>()[\]\\.,;:\s@"]+(\.[^<>()[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;
  return pattern.test(value) || "Invalid e-mail.";
};
