<template>
  <v-main>
    <v-container fluid fill-height>
      <v-row align="center" justify="center">
        <v-col xs="12" sm="8" md="4">
          <v-card elevation="12">
            <v-toolbar dark color="primary">
              <v-toolbar-title>{{ appName }} - Reset Password</v-toolbar-title>
            </v-toolbar>
            <v-card-text>
              <p class="subtitle-3">Enter your new password below</p>
              <v-form @keyup.enter="submit" v-model="valid" ref="form" @submit.prevent="" lazy-validation>
                <v-text-field
                  type="password"
                  ref="password"
                  label="Password"
                  :rules="password1Rules"
                  v-model="password1"
                ></v-text-field>
                <v-text-field
                  type="password"
                  label="Confirm Password"
                  :rules="password2Rules"
                  v-model="password2"
                ></v-text-field>
              </v-form>
            </v-card-text>
            <v-card-actions>
              <v-spacer></v-spacer>
              <v-btn @click="cancel">Cancel</v-btn>
              <v-btn @click="reset">Clear</v-btn>
              <v-btn @click="submit" :disabled="!valid">Save</v-btn>
            </v-card-actions>
          </v-card>
        </v-col>
      </v-row>
    </v-container>
  </v-main>
</template>

<script lang="ts">
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class UserProfileEdit extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  readonly password1Rules = [required];
  readonly password2Rules = [required, this.passwordIsEqual];

  passwordIsEqual(v) {
    return v === this.password1 || "Password should be the same";
  }

  appName = appName;
  valid = true;
  password1 = "";
  password2 = "";

  mounted() {
    this.checkToken();
  }

  reset() {
    this.password1 = "";
    this.password2 = "";
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.push("/");
  }

  checkToken() {
    const token = this.$router.currentRoute.query.token as string;
    if (!token) {
      this.mainContext.mutations.addNotification({
        content: "No token provided in the URL, start a new password recovery",
        color: "error",
      });
      this.$router.push("/password-recovery");
    } else {
      return token;
    }
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const token = this.checkToken();
      if (token) {
        await this.mainContext.actions.resetPassword({ token, password: this.password1 });
        this.$router.push("/");
      }
    }
  }
}
</script>
