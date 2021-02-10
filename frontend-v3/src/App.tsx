import React from "react";
import "./App.scss";
import { MainNavBar } from "./views/MainNavBar";
import { FocusStyleManager } from "@blueprintjs/core";
import { LoginForm } from "./views/LoginForm";
import { BrowserRouter } from "react-router-dom";
import { useAuthStore } from "./state/auth";

FocusStyleManager.onlyShowFocusOnTabs();

function App() {
  const loggedIn = useAuthStore(state => state.loggedIn);
  if (!loggedIn) {
    return <LoginForm />;
  }

  return (
    <BrowserRouter>
      <div className="bp3-dark">
        <MainNavBar />
        <main>Main</main>
      </div>
    </BrowserRouter>
  );
}

export default App;
