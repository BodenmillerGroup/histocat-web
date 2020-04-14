<template>
  <v-content>
    <v-container fluid fill-height>
      <v-row align="center" justify="center">
        <v-col xs="12" sm="8" md="4">
          <v-card elevation="12">
            <v-toolbar dark color="primary">
              <v-toolbar-title>Create {{ appName }} Account</v-toolbar-title>
              <v-spacer></v-spacer>
            </v-toolbar>
            <v-card-text>
              <v-form @keyup.enter="submit" v-model="valid" ref="form" @submit.prevent="" lazy-validation>
                <v-text-field
                  @keyup.enter="submit"
                  type="email"
                  v-model="email"
                  prepend-icon="mdi-account"
                  name="email"
                  label="Email"
                  :rules="emailRules"
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="text"
                  label="Full Name"
                  prepend-icon="mdi-account"
                  v-model="fullName"
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="password"
                  ref="password"
                  label="Password"
                  :rules="password1Rules"
                  v-model="password1"
                  prepend-icon="mdi-lock"
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="password"
                  label="Confirm Password"
                  :rules="password2Rules"
                  v-model="password2"
                  prepend-icon="mdi-lock"
                ></v-text-field>
              </v-form>
              <v-row>
                <v-col class="caption text-right py-0">
                  <router-link to="/login">Already have an account?</router-link>
                </v-col>
              </v-row>
            </v-card-text>
            <v-card-actions>
              <v-spacer></v-spacer>
              <v-btn @click.prevent="submit">Sign Up</v-btn>
            </v-card-actions>
          </v-card>
        </v-col>
      </v-row>
    </v-container>
  </v-content>
</template>

<script lang="ts">
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { userModule } from "@/modules/user";
import { email, required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class SignUp extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly userContext = userModule.context(this.$store);

  readonly emailRules = [required, email];
  readonly password1Rules = [required];
  readonly password2Rules = [required, this.passwordIsEqual];

  passwordIsEqual(v) {
    return v === this.password1 || "Password should be the same";
  }

  valid = true;
  email = "";
  fullName = "";
  password1 = "";
  password2 = "";
  appName = appName;

  reset() {
    this.email = "";
    this.fullName = "";
    this.password1 = "";
    this.password2 = "";
    (this.$refs.form as any).resetValidation();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const userExist = await this.checkUserExists();
      if (!userExist) {
        await this.userContext.actions.signUp({ email: this.email, password: this.password1 });
      }
    }
  }

  private async checkUserExists() {
    const userExist = await this.userContext.actions.checkUserExists(this.email);
    if (userExist) {
      this.mainContext.mutations.addNotification({
        content: "User with this email already exists",
        color: "warning",
      });
    }
    return userExist;
  }
}
</script>

<style></style>
