import React, { useEffect } from "react";
import { MainNavBar } from "views/MainNavBar";
import { FocusStyleManager } from "@blueprintjs/core";
import { LoginForm } from "views/auth/LoginForm";
import { BrowserRouter, Switch, Route } from "react-router-dom";
import { useAuthStore } from "modules/auth";
import history from "utils/history";
import { SignupForm } from "./views/auth/SignupForm";
import { PasswordRecoveryForm } from "./views/auth/PasswordRecoveryForm";
import { ResetPasswordForm } from "./views/auth/ResetPasswordForm";
import { ProfileView } from "./views/profile/ProfileView";
import "./App.scss";
import { ProfileEditView } from "./views/profile/ProfileEditView";
import { ProfilePasswordView } from "./views/profile/ProfilePasswordView";
import { AdminView } from "./views/admin/AdminView";
import { UsersView } from "./views/admin/users/UsersView";

FocusStyleManager.onlyShowFocusOnTabs();

function App() {
  const { loggedIn, checkLoggedIn } = useAuthStore();

  useEffect(() => {
    checkLoggedIn();
  }, [checkLoggedIn]);

  if (!loggedIn) {
    switch (history.location.pathname) {
      case "/signup":
        return <SignupForm />;
      case "/password-recovery":
        return <PasswordRecoveryForm />;
      case "/reset-password":
        return <ResetPasswordForm />;
      default:
        return <LoginForm />;
    }
  }

  return (
    <BrowserRouter>
      <div className="bp3-dark">
        <MainNavBar />
        <main>
          <Switch>
            <Route exact={true} path="/profile">
              <ProfileView />
            </Route>
            <Route path="/profile/edit">
              <ProfileEditView />
            </Route>
            <Route path="/profile/password">
              <ProfilePasswordView />
            </Route>

            <Route path="/admin">
              <AdminView/>
            </Route>
          </Switch>
        </main>
      </div>
    </BrowserRouter>
  );
}

export default App;
