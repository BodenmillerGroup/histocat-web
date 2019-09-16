<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Edit User</div>
      </v-card-title>
      <v-card-text>
        <template>
          <div class="my-4">
            <div class="subtitle-1 primary--text text--lighten-2">Username</div>
            <div class="title primary--text text--darken-2" v-if="user">{{ user.email }}</div>
            <div class="title primary--text text--darken-2" v-else>-----</div>
          </div>
          <v-form v-model="valid" ref="form">
            <v-text-field label="Full Name" v-model="fullName"></v-text-field>
            <v-text-field label="E-mail" type="email" v-model="email" :rules="emailRules"></v-text-field>
            <v-checkbox label="Is Superuser" v-model="isSuperuser"></v-checkbox>
            <v-checkbox label="Is Active" v-model="isActive"></v-checkbox>
            <v-row align="center">
              <v-col class="shrink">
                <v-checkbox v-model="setPassword"></v-checkbox>
              </v-col>
              <v-col>
                <v-text-field
                  :disabled="!setPassword"
                  type="password"
                  ref="password"
                  label="Set Password"
                  :rules="password1Rules"
                  v-model="password1"
                >
                </v-text-field>
                <v-text-field
                  v-show="setPassword"
                  type="password"
                  label="Confirm Password"
                  :rules="password2Rules"
                  v-model="password2"
                >
                </v-text-field>
              </v-col>
            </v-row>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { userModule } from "@/modules/user";
import { IUserProfileUpdate } from "@/modules/user/models";
import { email, required } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class EditUser extends Vue {
  readonly userContext = userModule.context(this.$store);

  readonly emailRules = [required, email];
  readonly password1Rules = [this.passwordRequired];
  readonly password2Rules = [this.passwordRequired, this.passwordIsEqual];

  passwordRequired(v) {
    return this.setPassword && !v ? "Required" : true;
  }

  passwordIsEqual(v) {
    return this.setPassword && v !== this.password1 ? "Password should be the same" : true;
  }

  valid = true;
  fullName: string = "";
  email: string = "";
  isActive: boolean = true;
  isSuperuser: boolean = false;
  setPassword = false;
  password1: string = "";
  password2: string = "";

  async mounted() {
    await this.userContext.actions.getUsers();
    this.reset();
  }

  reset() {
    this.setPassword = false;
    this.password1 = "";
    this.password2 = "";
    if (this.user) {
      this.fullName = this.user.full_name;
      this.email = this.user.email;
      this.isActive = this.user.is_active;
      this.isSuperuser = this.user.is_superuser;
    }
    if (this.$refs.form) {
      (this.$refs.form as any).resetValidation();
    }
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const data: IUserProfileUpdate = {};
      if (this.fullName) {
        data.full_name = this.fullName;
      }
      if (this.email) {
        data.email = this.email;
      }
      data.is_active = this.isActive;
      data.is_superuser = this.isSuperuser;
      if (this.setPassword) {
        data.password = this.password1;
      }
      await this.userContext.actions.updateUser({ id: this.user!.id, user: data });
      this.$router.push("/main/admin/users");
    }
  }

  get user() {
    return this.userContext.getters.getUser(+this.$router.currentRoute.params.id);
  }
}
</script>
